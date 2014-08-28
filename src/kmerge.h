#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <pthread.h>
#include "H5Cpp.h"
#include "fq.h"
#include "hdf5file.h"
#include "SpookyV2.h"
#include <dlib/logger.h> 
#include "klib/kseq.h"
#include <zlib.h>
#include "blosc_filter.h"
#include <stx/btree_map>
#include "cpp-btree/btree_map.h"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

KSEQ_INIT(gzFile, gzread)

#define MAX_UINT_VAL 4294967295 //2^32-1

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
  std::string compound_name;
  std::string hash_dataset_name;
  std::string counts_dataset_name;
  uint num_threads;
  uint priority;
  std::string lock_filename;
  std::vector<std::string> serialized_files;
} ;


class KMerge {
 private:
  HDF5 *hdf5_file;
  HashEnumType hash_type;
  std::string dir;
  static pthread_mutex_t mutex;
  dlib::logger dlog;
  std::string filename;

 public:
  static const uint CHUNK_ROW_SIZE = 2097152; // calculated this based off of pytables optimization of parameter

  KMerge(const std::string&, const std::string&, const std::string&);
  ~KMerge();
  static std::string rev_comp(const std::string&);
  void build(param_struct&);
  bool count_hashed_kmers(param_struct&, bool);
  bool add_dataset(const uint, const uint*, param_struct&);
  bool add_dataset(const std::string, uint, const uint*, HDF5*);
  bool add_taxonomy(const std::string&, const std::string&);
  static uint hash_kmer(const std::string&, const HashEnumType);
  uint hash_kmer(const std::string&);
  bool add_hash_and_count(std::vector<uint>&, std::vector<uint>&, uint, uint);
  bool add_hash_and_count(std::map<uint, uint>&, uint, uint);
  bool add_hash(btree::btree_map<uint, uint>&, uint);
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


  class CountAndHashSeqFastx {
  public:

    param_struct& params;
    bool print_status;

    void operator() (long i) const;

  CountAndHashSeqFastx( param_struct& params_, bool print_status_) : params(params_), print_status(print_status_) {}

    ~CountAndHashSeqFastx() {}
    };

};

