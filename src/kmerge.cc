#include "kmerge.h"
#include "lookup3.c"
#include "SpookyV2.h"
#include "MurmurHash3.h"
#include "city.h"
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <fstream>

using namespace seqan;
using namespace std;

pthread_mutex_t KMerge::mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t KMerge::mem_mutex = PTHREAD_MUTEX_INITIALIZER;

KMerge::KMerge (const std::string& filename, const std::string& hash_func, const std::string& dir) {
  this->hdf5_file = new HDF5(filename, false);
  this->dir = dir;

  if (hash_func == "lookup3") {
    this->hash_type = LOOKUP3;
  } else if (hash_func == "spooky") {
    this->hash_type = SPOOKY;
  } else if (hash_func == "murmur") {
    this->hash_type = MURMUR;
  } else if (hash_func == "city") {
    this->hash_type = CITY;
  } else {
    throw "Invalid hash function provided";
  }

}

KMerge::~KMerge() {
  delete this->hdf5_file;
}

uint KMerge::hash_kmer(const std::string& kmer, const HashEnumType hash_type) {
  string rc_kmer(kmer.c_str()), combine;
  reverseComplement(rc_kmer);
  if (kmer < rc_kmer) {
    combine = kmer + rc_kmer;
  } else {
    combine = rc_kmer + kmer;
  }

  switch (hash_type) {
  case LOOKUP3:
    return(hashlittle(combine.c_str(), combine.length(), 0));
  case SPOOKY:
    return(SpookyHash::Hash32(combine.c_str(), combine.length(), 0));
  case MURMUR:
    uint out;
    MurmurHash3_x86_32  ( combine.c_str(), combine.length(), 0, &out );
    return out;
  case CITY:
    return(CityHash32(combine.c_str(), combine.length()));
  default:
    throw "Invalid hash function provided";
  }
}

uint KMerge::hash_kmer(const std::string& kmer) {
  return KMerge::hash_kmer(kmer, this->hash_type);
}

bool KMerge::add_hash_and_count(std::vector<uint>& hashes, std::vector<uint>& counts, uint kmer_hash_val, uint kmer_count) {
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

bool KMerge::add_hash_and_count(std::map<uint, uint>& hashed_counts, uint kmer_hash_val, uint kmer_count) {
  if (hashed_counts.find(kmer_hash_val) == hashed_counts.end()) {
    hashed_counts[kmer_hash_val] = kmer_count;
  } else {
    hashed_counts[kmer_hash_val] += kmer_count;
  }
  return true;
}

bool KMerge::add_hash(std::map<uint, uint>& hashed_counts, uint kmer_hash_val) {
  if (hashed_counts.find(kmer_hash_val) == hashed_counts.end()) {
    hashed_counts[kmer_hash_val] = 1; 
  } else {   
    hashed_counts[kmer_hash_val]++; 
  }
  return true;
}

bool KMerge::count_hashed_kmers(std::string& filename, uint k, std::map<uint, uint>& hashed_counts) {

  seqan::SequenceStream seq_io(filename.c_str());

  seqan::CharString id;
  seqan::Dna5String seq;
  std::stringstream kmer;

  while (!atEnd(seq_io)) {
    if (readRecord(id, seq, seq_io) != 0) continue; 
    if (length(seq) < k) continue;
    for (int i = 0; i < length(seq) - k + 1; i++) {
      seqan::Infix<seqan::Dna5String>::Type sub_seq(seq, i, i+k);

      kmer.str("");
      kmer << sub_seq;

      if(kmer.str().find("N") != std::string::npos) { // skip kmers containing Ns
        continue;
      }
      if(!this->add_hash(hashed_counts, this->hash_kmer(kmer.str()))) { 
	throw "Unable to add hash and count";
	return false;
      }
    }
  }

  return true;
}


bool KMerge::add_dataset(const std::string dataset_path, const uint data_size, const uint* data, HDF5* h5_file=NULL) {
  std::vector<uint64_t> dims;

  if (h5_file == NULL) {
    h5_file = this->hdf5_file;
  }

  dims.push_back(data_size);
  if (!(h5_file->createDataset(dataset_path, dims, FQ::FQT_UINT))) {
    throw "Unable to create dataset";
    return false;
  }
  if (!(h5_file->setData(dataset_path, data))) {
    throw "Unable to set data for dataset";
    return false;
  }

  return true;
}

uint* KMerge::get_dataset(const std::string h5_filename, const uint data_size, const std::string dataset_path) {
  HDF5 *h5 = new HDF5(h5_filename);

  uint* data = new uint[data_size];

  if(!(h5->getData(dataset_path, data))) {
    cerr << "Cannot access data" << endl;
    exit(EXIT_FAILURE);
  }
  
  delete h5;
  
  return data;
}



bool KMerge::add_taxonomy(const std::string& group) {
  std::stringstream path, in_file_ss, path_root, error;
  std::vector<std::string> lines;
  std::vector<uint64_t> dims;
  std::string line;


  path_root << group << "/taxonomy";
  in_file_ss << this->dir << path_root.str() << ".txt";
  ifstream in_file(in_file_ss.str().c_str());


  dims.push_back(0);

  while (std::getline(in_file, line)) {
    std::istringstream tokenizer(line);
    path.str("");
    path << path_root.str();
    while (!tokenizer.eof()) {
      std::string token;
      getline(tokenizer, token, '\t');
      path << "/" << token;
    }
    if (!(this->hdf5_file->createDataset(path.str(), dims, FQ::FQT_BYTE))) {
      error << "Unable to add classification:" << path.str();
      throw error.str();
      return false;
    }
  }

  return true;
}

bool KMerge::sort_kmer_hashes_and_counts(std::vector<uint>& hashes, std::vector<uint>& counts) {
  std::map<uint, uint> sorter;

  std::transform( hashes.begin(), hashes.end(), counts.begin(),
		  std::inserter(sorter, sorter.end() ), std::make_pair<uint,int> );
  
  //clear hashes and counts and associated memory
  std::vector<uint>().swap( hashes );
  std::vector<uint>().swap( counts );

  for (std::map<uint, uint>::const_iterator iter = sorter.begin(); iter != sorter.end(); iter++) {
    hashes.push_back(iter->first);
    counts.push_back(iter->second);
  }
  return true;
}

void KMerge::BuilderTask::execute() {
  stringstream file_name, file_loc, error;
  std::map<uint, uint> hashed_counts;

  cout << "Working on " << params.group_name << endl;
  for (uint k = params.k_val_start; k <= params.k_val_end; k=k+2) {
    //if k = optimal k, also get raw k-mer sequence
    if(!(params.kmerge->count_hashed_kmers(params.seq_filename, k, hashed_counts))) {
      error << "Unable to parse " << params.seq_filename << " for k=" << k << std::endl;
      throw error.str();
    } else {
      cout << "Finished parsing: " << params.seq_filename << " for k=" << k << std::endl;
      cout << "Hashes vector size now: " << hashed_counts.size() << std::endl;  
    }
    error.str("");
  }

  uint hash_count = hashed_counts.size();

  std::vector<uint> hashes(hash_count), counts(hash_count);

  uint i = 0;
  for (map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    hashes[i] = iter->first;
    counts[i] = iter->second;
    i++;
  }

  // remove all elements from map as they are no longer needed
  std::map<uint, uint>().swap( hashed_counts );

  HDF5 *out_h5 = new HDF5(params.tmp_h5_filename, false);

  if(!(params.kmerge->add_dataset("/kmer_hash", hash_count, &hashes[0], out_h5))) {
    error << "Unable to add hashes for " << params.group_name << " to temp file." << endl;
    throw error.str();
  }

  //clear memory
  std::vector<uint>().swap(hashes);

  if(!(params.kmerge->add_dataset("/count", hash_count, &counts[0], out_h5))) {
    error << "Unable to add counts for " << params.group_name << " to temp file." << endl;
    throw error.str();
  }

  //clear memory
  std::vector<uint>().swap(counts);

  delete out_h5;

  pthread_mutex_lock( &KMerge::mutex );

  uint* hashes_arr = params.kmerge->get_dataset(params.tmp_h5_filename, hash_count, "/kmer_hash");

  if(!(params.kmerge->add_dataset(params.hash_dataset_name, hash_count, hashes_arr, NULL))) {
    error << "Unable to add hashes for " << params.group_name << endl;
    throw error.str();
  }

  delete [] hashes_arr;

  uint* counts_arr = params.kmerge->get_dataset(params.tmp_h5_filename, hash_count, "/count");

  if(!(params.kmerge->add_dataset(params.counts_dataset_name, hash_count, counts_arr, NULL))) {
    error << "Unable to add counts for " << params.group_name << endl;
    throw error.str();
  }

  delete [] counts_arr;

  if(!(params.kmerge->add_taxonomy(params.group_name))) {
    error << "Unable to add classifications for " << params.group_name << endl;
    throw error.str();
  }

  pthread_mutex_unlock( &KMerge::mutex );

  cout << hash_count << " k-mer hashes for " << params.group_name << endl;

  cout << "Done (" << params.group_name  << ")" << endl;

  remove(params.tmp_h5_filename.c_str());

  return;

}
