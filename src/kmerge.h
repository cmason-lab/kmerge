#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <pthread.h>
#include "SpookyV2.h"
#include <dlib/logger.h> 
#include "klib/kseq.h"
#include <zlib.h>
#include "cpp-btree/btree_map.h"
#include <ulib/hash_chain.h>
#include "leveldb/db.h"

KSEQ_INIT(gzFile, gzread)

#define MAX_UINT_VAL 4294967295 //2^32-1

using namespace std;

enum HashEnumType { LOOKUP3, SPOOKY, MURMUR, CITY };

typedef unsigned int uint;

class KMerge;

struct param_struct {
  KMerge * kmerge;
  std::string db_filename;
  uint k_val_start;
  uint k_val_end;
  std::string seq_filename;
  std::string group_name;
  std::string compound_name;
  uint num_threads;
  uint priority;
  std::string lock_filename;
  leveldb::DB* db;
} ;


class KMerge {
 private:
  HashEnumType hash_type;
  std::string dir;
  dlib::logger dlog;
  std::string filename;
  static pthread_mutex_t db_mutex;

 public:
  static const uint CHUNK_ROW_SIZE = 2097152; // calculated this based off of pytables optimization of parameter
  static const uint INIT_MAP_CAPACITY = 100000000; //used to initialize chain_hash_map

  KMerge(const std::string&, const std::string&, const std::string&);
  ~KMerge();
  static std::string rev_comp(const std::string&);
  static std::vector<uint> compress(const std::vector<uint>&);
  static std::vector<uint> uncompress(const std::vector<uint>&, uint);
  void build(param_struct&);
  bool count_hashed_kmers(param_struct&, ulib::chain_hash_map<uint, uint>&, bool);
  bool add_dataset_size(uint, const std::string&, leveldb::DB*);
  bool add_dataset(const std::vector<uint>&, const std::string&, leveldb::DB*);
  bool add_dataset(const uint, const uint*, param_struct&);
  bool add_taxonomy(const std::string&, leveldb::DB*);
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

  class CountAndHashSeq {
  public:

    param_struct& params;
    ulib::chain_hash_map<uint, uint>& hashed_counts;
    bool print_status;

    void operator() (long i) const;

  CountAndHashSeq( param_struct& params_, ulib::chain_hash_map<uint, uint>& hashed_counts_, bool print_status_) : params(params_), hashed_counts(hashed_counts_), print_status(print_status_) {}

    ~CountAndHashSeq() {}
    };

};

