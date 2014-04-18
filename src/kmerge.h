#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <pthread.h>
#include "fq.h"
#include "hdf5file.h"


#define MAX_UINT_VAL 4294967295 //2^32-1

using namespace std;

typedef unsigned int uint;

class KMerge;

struct param_struct {
  KMerge * kmerge;
  std::string hdf5_filename;
  uint k_val_start;
  uint k_val_end;
  std::string seq_filename;
  std::string group_name;
  std::string hash_dataset_name;
  std::string counts_dataset_name;
} ;



class KMerge {
 private:
  std::string filename;
  HDF5 *hdf5_file;
  static pthread_mutex_t mutex;
  std::string hash_function;

 public:
  static const uint GZIP_BEST_COMPRESSION = 9;


  KMerge(const std::string&, const std::string&);
  ~KMerge();
  static uint hash_kmer(const std::string&, const std::string&);
  uint hash_kmer(const std::string&);
  bool count_hashed_kmers(std::string&, uint, std::map<uint, uint>&);
  static bool add_hash_and_count(std::vector<uint>&, std::vector<uint>&, uint, uint);
  static bool add_hash_and_count(std::map<uint, uint>&, uint, uint);
  bool add_dataset(const std::string, uint, const uint*);
  static bool sort_kmer_hashes_and_counts(std::vector<uint>&, std::vector<uint>&);

  class BuilderTask {
    public:
      param_struct params;

      void execute();
    
      BuilderTask(const param_struct params) {
	this->params = params;
      }

      ~BuilderTask() {
      }
  };
};
