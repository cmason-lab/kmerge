#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <pthread.h>
#include "fq.h"
#include "hdf5file.h"


#define MAX_UINT_VAL 4294967295 //2^32-1
#define THROTTLE_KMER_LENGTH 19 //lock k-mer counting above this value to throttle memory allocation for longer k-mers

using namespace std;

enum HashEnumType { LOOKUP3, SPOOKY, MURMUR, CITY };

typedef unsigned int uint;

class KMerge;

struct param_struct {
  KMerge * kmerge;
  std::string hdf5_filename;
  uint k_val_start;
  uint k_val_end;
  std::string seq_filename;
  std::string tmp_hashes_filename;
  std::string tmp_counts_filename;
  std::string group_name;
  std::string hash_dataset_name;
  std::string counts_dataset_name;
} ;


class KMerge {
 private:
  HDF5 *hdf5_file;
  HashEnumType hash_type;
  std::string dir;
  static pthread_mutex_t mutex;

 public:
  KMerge(const std::string&, const std::string&, const std::string&);
  ~KMerge();
  bool count_hashed_kmers(std::string&, uint, std::map<uint, uint>&);
  bool add_dataset(const std::string, uint, const uint*, HDF5*);
  bool add_taxonomy(const std::string&);
  static uint hash_kmer(const std::string&, const HashEnumType);
  uint hash_kmer(const std::string&);
  bool add_hash_and_count(std::vector<uint>&, std::vector<uint>&, uint, uint);
  bool add_hash_and_count(std::map<uint, uint>&, uint, uint);
  bool add_hash(std::map<uint, uint>&, uint);
  bool sort_kmer_hashes_and_counts(std::vector<uint>&, std::vector<uint>&);

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

