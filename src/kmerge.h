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
#include "unqlite.h"
#include "fastpfor/memutil.h"
#include <mutex>     
#include <condition_variable>


KSEQ_INIT(gzFile, gzread)

#define MAX_UINT_VAL 4294967295 //2^32-1
#define BYTES_IN_GB 1073741824.0
#define BYTE_ALIGNED_SIZE 16
#define MIN_MEM 0 //GB
#define POLL_INTERVAL 5 //seconds

using namespace std;


enum HashEnumType { LOOKUP3, SPOOKY, MURMUR, CITY };

typedef unsigned int uint;

class KMerge;

struct param_struct {
  KMerge * kmerge;
  std::string db_filename;
  std::string dump_filename;
  uint k_val_start;
  uint k_val_end;
  std::string seq_filename;
  std::string group_name;
  std::string compound_name;
  uint num_threads;
  uint priority;
  std::string lock_filename;
  std::mutex* dump_mtx;
  std::condition_variable* cv;
  bool ready;
  bool finished_hashing;
  unqlite* db;
  ulib::chain_hash_map<uint, uint>* hashed_counts;
  bool is_ref;
  std::vector<char> writing;
  bool polling_done;
} ;


class KMerge {
 private:
  HashEnumType hash_type;
  std::string dir;
  dlib::logger dlog;
  std::string filename;
  static pthread_mutex_t db_mutex;
  double max_gb;

 public:
  static const uint CHUNK_ROW_SIZE = 2097152; // calculated this based off of pytables optimization of parameter
  static const uint INIT_MAP_CAPACITY = 100000000; //used to initialize chain_hash_map
  static const uint PARTITION_SIZE = 500000000;

  KMerge(const std::string&, const std::string&, const std::string&, double max_gb=MIN_MEM);
  ~KMerge();
  static std::string rev_comp(const std::string&);
  static std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > compress(const std::vector<uint>&);
  static std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > uncompress(const std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >&, uint); 
  static double memory_used(void);
  static void dump_hashes(ulib::chain_hash_map<uint, uint>&, std::string&);
  static void load_hashes(btree::btree_map<uint, uint>&, std::string&);
  void poll_memory(param_struct&);
  void build(param_struct&);
  bool count_hashed_kmers(param_struct&, ulib::chain_hash_map<uint, uint>&, bool);
  bool add_property(uint, const std::string&, unqlite*);
  bool add_dataset(const std::vector<uint>&, const std::string&, unqlite*, bool);
  bool add_dataset(const uint, const uint*, param_struct&);
  bool add_taxonomy(const std::string&, unqlite*);
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
      void check_memory();
    
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

