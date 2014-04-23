#include "kmerge.h"
#include "lookup3.c"
#include "SpookyV2.h"
#include "MurmurHash3.h"
#include "city.h"
#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>
#include <fstream>
#include <libGkArrays/gkArrays.h>

using namespace seqan;
using namespace std;

pthread_mutex_t KMerge::mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t KMerge::mem_mutex = PTHREAD_MUTEX_INITIALIZER;

KMerge::KMerge (const std::string& filename, const std::string& hash_func) {
  this->hdf5_file = new HDF5(filename, false);
  
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

uint KMerge::hash_kmer(const std::string& kmer) {
  string rc_kmer(kmer.c_str()), combine;
  reverseComplement(rc_kmer);
  if (kmer < rc_kmer) {
    combine = kmer + rc_kmer;
  } else {
    combine = rc_kmer + kmer;
  }

  switch (this->hash_type) {
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

bool KMerge::count_hashed_kmers(std::string& filename, uint k, std::map<uint, uint>& hashed_counts /*, get_raw=false */) {
  std::set<std::string> kmers;
  // Building the index
  gkarrays::gkArrays *genome = new gkarrays::gkArrays(&filename[0], k, true /*use bitvector to conserve space*/, 
						      0 /*not specifying read length*/, 
						      true /* stranded counting (rev comp counted seperately)*/);
  gkarrays::readIterator *read_iter = genome->getReads()->begin();
  uint i = 0, num_iterations = 0;
  bool all_kmers_found = false;
  while (!read_iter->isFinished() && !all_kmers_found) {
    uint *support = genome->getSupport(i);
    for(uint j = 0; j < genome->getSupportLength(i); j++) {
      num_iterations++;
      char *kmer = genome->getTagFactor(i, j, k);
      if (kmers.find(kmer) == kmers.end()) {
	kmers.insert(kmer);
	if(!this->add_hash_and_count(hashed_counts, this->hash_kmer(kmer), support[j])) {
	  throw "Unable to add hash and count";
	  return false;
	}
      }
      if (genome->getGkCFALength() == kmers.size()) { //we've seen all kmers in the genome
        all_kmers_found = true;
	delete [] kmer;
        break;
      }
      delete [] kmer;
    }
    delete [] support;
    ++(*read_iter);
    i++;
  }

  cout << "Found " << kmers.size() << " distinct " << k << "-mers" << std::endl;

  delete read_iter;
  delete genome;

  return true;
}

bool KMerge::add_dataset(const std::string dataset_path, uint data_size, const uint* data) {
  std::vector<uint64_t> dims;


  dims.push_back(data_size);
  if (!(this->hdf5_file->createDataset(dataset_path, dims, FQ::FQT_UINT))) {
    throw "Unable to create dataset";
    return false;
  }
  if (!(this->hdf5_file->setData(dataset_path, data))) {
    throw "Unable to set data for dataset";
    return false;
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
  std::vector<uint> hashes;
  std::vector<uint> counts;

  cout << "Working on " << params.group_name << endl;
  for (uint k = params.k_val_start; k <= params.k_val_end; k=k+2) {
    if (k > THROTTLE_KMER_LENGTH) pthread_mutex_lock( &KMerge::mem_mutex ); //throttle memory for longer k-mers
    //if k = optimal k, also get raw k-mer sequence
    if(!(params.kmerge->count_hashed_kmers(params.seq_filename, k, hashed_counts))) {
      error << "Unable to parse " << params.seq_filename << " for k=" << k << std::endl;
      throw error.str();
    } else {
      cout << "Finished parsing: " << params.seq_filename << " for k=" << k << std::endl;
      cout << "Hashes vector size now: " << hashed_counts.size() << std::endl;  
    }
    if (k > THROTTLE_KMER_LENGTH) pthread_mutex_unlock( &KMerge::mem_mutex );
    error.str("");
  }
  for (map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    hashes.push_back(iter->first);
    counts.push_back(iter->second);
  }
  // remove all elements from map as they are no longer needed
  std::map<uint, uint>().swap( hashed_counts );

  pthread_mutex_lock( &KMerge::mutex );

  if(!(params.kmerge->add_dataset(params.hash_dataset_name, hashes.size(), &hashes[0]))) {
    error << "Unable to add hashes for " << params.group_name << endl;
    throw error.str();
  }
  if(!(params.kmerge->add_dataset(params.counts_dataset_name, counts.size(), &counts[0]))) {
    error << "Unable to add counts for " << params.group_name << endl;
    throw error.str();
  }
  pthread_mutex_unlock( &KMerge::mutex );

  cout << hashes.size() << " k-mer hashes for " << params.group_name << endl;
  // clear memory
  std::vector<uint>().swap(hashes);
  std::vector<uint>().swap(counts);
  cout << "Done (" << params.group_name  << ")" << endl;
  return;

}
