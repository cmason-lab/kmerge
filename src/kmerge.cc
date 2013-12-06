#include "kmerge.h"
#include "lookup3.c"
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>
#include <fstream>

using namespace seqan;
using namespace std;

#define ROW_SIZE 0
#define COL_SIZE 1

pthread_mutex_t KMerge::mutex = PTHREAD_MUTEX_INITIALIZER;

KMerge::KMerge (const H5std_string& file_name) {
  this->file_name = file_name;
}

uint KMerge::hashKmer(const std::string& kMer) {
  return((uint) hashlittle(kMer.c_str(), kMer.length(), 0));
}

bool KMerge::addHashAndCount(std::vector<uint>& hashes, std::vector<uint>& counts, uint kmer_hash_val, uint kmer_count) {
  std::vector<uint>::iterator iter;
  uint pos = 0;                                                                                                                    
  for (iter=hashes.begin(); iter!=hashes.end(); iter++) {                                                                                                                                             
    if (*iter==kmer_hash_val) {                                                                                                                                                              
      break;                                                                                                                                                                                                        
    }                                                                                                                                                                                                              
    pos++;                                                                                           
  }
  if (iter == hashes.end()) {
    hashes.push_back(kmer_hash_val);                                                                                                                   
    counts.push_back(kmer_count);
  } else {
    counts[pos] += kmer_count;
  }
  return true;
}

bool KMerge::parseKmerCountsFile(const std::string& counts_file_name, std::vector<uint>& hashes, std::vector<uint>& counts) {
  Stream<GZFile> gzStream;
  std::string kMer;
  std::string sCount;

  if (!open(gzStream, counts_file_name.c_str(), "r")) {
    cerr << "ERROR: Could not open file " << counts_file_name << '\n';        
    return false;
  }

  RecordReader<Stream<GZFile>, SinglePass<> > reader(gzStream);
    
  while (!atEnd(reader)) {
    clear(kMer);
    if (readUntilWhitespace(kMer, reader) != 0) {
      cerr << "Problem with file" << endl;
      return false;
    }
    goNext(reader);
    sCount.clear();
    if (readUntilWhitespace(sCount, reader) != 0) {
      cerr << "Problem with file" << endl;
      return false;
    }
    goNext(reader);
    this->addHashAndCount(hashes, counts, KMerge::hashKmer(kMer), atoi(sCount.c_str()));
  }
  close(gzStream);

  return true;

}

bool KMerge::addHashAndCount(std::map<uint, uint>& hashed_counts, uint kmer_hash_val, uint kmer_count) {
  std::map<uint, uint>::iterator iter;

  iter = hashed_counts.find(kmer_hash_val);
  
  if (iter == hashed_counts.end()) {
    hashed_counts[kmer_hash_val] = kmer_count;
  } else {
    uint prev_count = iter->second;
    hashed_counts[kmer_hash_val] = prev_count + kmer_count;
  }
  return true;
}

bool KMerge::parseKmerCountsFile(const std::string& counts_file_name, std::map<uint, uint>& hashed_counts) {
  Stream<GZFile> gzStream;
  std::string kMer;
  std::string sCount;

  if (!open(gzStream, counts_file_name.c_str(), "r")) {
    cerr << "ERROR: Could not open file " << counts_file_name << '\n';        
    return false;
  }

  RecordReader<Stream<GZFile>, SinglePass<> > reader(gzStream);
    
  while (!atEnd(reader)) {
    clear(kMer);
    if (readUntilWhitespace(kMer, reader) != 0) {
      cerr << "Problem with file" << endl;
      return false;
    }
    goNext(reader);
    sCount.clear();
    if (readUntilWhitespace(sCount, reader) != 0) {
      cerr << "Problem with file" << endl;
      return false;
    }
    goNext(reader);
    this->addHashAndCount(hashed_counts, KMerge::hashKmer(kMer), atoi(sCount.c_str()));
  }
  close(gzStream);

  return true;

}

bool KMerge::addDatasetToHDF5File(const H5std_string& group_name, const H5std_string& ds_name, const hsize_t data_size, const uint* data, const bool create_group = false) {
  ifstream ifile(this->file_name.c_str());
  if (ifile) {
    this->file = new H5File( this->file_name, H5F_ACC_RDWR );
  } else {
    this->file = new H5File( this->file_name, H5F_ACC_TRUNC );
  }

  /*                              
   * Create property list for a dataset and set up fill values. 
   */
  uint rank = 2;
  uint fillvalue = 0;   /* Fill value for the dataset */
  DSetCreatPropList plist;
  hsize_t max_dims[2] = {H5S_UNLIMITED, 1};
  Group* group = NULL;
  hsize_t chunk_dims[2] = {KMerge::CHUNK_ROW_SIZE, 1};
  plist.setChunk(rank, chunk_dims);
  plist.setFillValue(PredType::NATIVE_UINT, &fillvalue);
  plist.setDeflate( KMerge::GZIP_BEST_COMPRESSION );

  /* 
   * Create dataspace for the dataset                   
   */
  const hsize_t data_dims[2] = {data_size, 1}; 
                                                                                                     
  DataSpace file_space( rank, data_dims, max_dims);

  /*
   * Add group
   */

  if (create_group) {
    group = new Group( this->file->createGroup( group_name ));
  }


  /* 
   * Create dataset and write it into the file.                        
   */

  DataSpace mem_space( rank, data_dims, max_dims );

  DataSet* dataset = new DataSet( file->createDataSet(ds_name, PredType::NATIVE_UINT, file_space, plist));
  dataset->write( data, PredType::NATIVE_UINT, mem_space, file_space );
  
  /*                                                                                                       
   * Close the dataset and free memory.
   */
  if (create_group) {
    group->close();
    delete group;
  }
  dataset->close();
  delete dataset;

  this->file->close();
  delete this->file;

  return true;
}

void KMerge::parseAndWriteInThread(void* arg) {
  param_struct * params = (param_struct*) arg;
  stringstream file_name, file_loc;
  map<uint, uint> hashed_counts;
  vector<uint> hashes;
  vector<uint> counts;
 
  cout << "Working on " << params->group_name << endl;

  for (uint k = params->k_val_start; k <= params->k_val_end; k=k+2) {
    file_name.str("");
    file_name << "k" << k << ".counts.gz";
    file_loc.str("");
    file_loc << params->group_name << "/" << file_name.str();
    
    if (!(params->kmerge->parseKmerCountsFile(file_loc.str(), hashed_counts))) {
      cerr << "Unable to parse: " << file_loc.str() << endl;
    } else {
      cout << "Finished parsing: " << file_loc.str() << endl;
      cout << "Hashes vector size now: " << hashed_counts.size() << endl;
    }
  }
  for (map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    hashes.push_back(iter->first);
    counts.push_back(iter->second);
  }
  // remove all elements from map as they are no longer needed
  hashed_counts.clear();
  pthread_mutex_lock( &KMerge::mutex );
  if(!(params->kmerge->addDatasetToHDF5File(params->group_name, params->hash_dataset_name, hashes.size(), &hashes[0], true))) {
    cerr << "Unable to add hashes for " << params->group_name << endl;
  }
  if(!(params->kmerge->addDatasetToHDF5File(params->group_name, params->count_dataset_name, counts.size(), &counts[0], false))) {
    cerr << "Unable to add counts for " << params->group_name << endl;
  }
  pthread_mutex_unlock( &KMerge::mutex );
  cout << hashes.size() << " k-mer hashes for " << params->group_name << endl;
  hashes.clear();
  counts.clear();
  cout << "Done (" << params->group_name  << ")" << endl;

  return;
}

KMerge::~KMerge () {
}
