#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "SpookyV2.h"
#include <dlib/logger.h> 
#include <dlib/compress_stream.h>
#include "klib/kseq.h"
#include <zlib.h>
#include "cpp-btree/btree_map.h"
#include <mutex>

KSEQ_INIT(gzFile, gzread)

#define MAX_UINT_VAL 4294967295 //2^32-1

using namespace std;

typedef dlib::compress_stream::kernel_3b cs;

enum HashEnumType { LOOKUP3, SPOOKY, MURMUR, CITY};
enum TruncType { BITS_8, BITS_16, BITS_32 };

typedef unsigned int uint;

class KMerge;

struct param_struct {
  KMerge * kmerge;
  uint k_val_start;
  uint k_val_end;
  std::string seq_filename;
  std::string group_name;
  uint num_threads;
  std::string hashes_filename;
  std::string counts_filename;
  bool is_ref;
};


class KMerge {
 private:
  HashEnumType hash_type;
  TruncType trunc_type;
  std::string dir;
  std::string out_dir;
  dlib::logger dlog;
  std::string filename;

 public:

  KMerge(const std::string&, const std::string&, const std::string&, const uint);
  ~KMerge();
  static std::string rev_comp(const std::string&);
  static uint hash_kmer(const std::string&, const HashEnumType, const TruncType);
  static std::string get_seq_base_id(const std::string&, const std::string&);
  bool count_hashed_kmers(param_struct&, btree::btree_map<uint, uint>&, bool, bool);
  bool hash_seq(std::string&, uint, btree::btree_map<uint, uint>&, std::mutex&);
  bool hash_seq(const std::vector<std::string>&, uint, btree::btree_map<std::string, btree::btree_map<uint, uint> >&, const std::string&, std::mutex&);
  bool add_dataset(const std::vector<uint>&, const std::string&);
  bool add_taxonomy(const std::string&);
  uint hash_kmer(const std::string&);


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

  class HashSeq {
  public:
    param_struct& params;
    btree::btree_map<uint, uint>& hashed_counts;
    const std::vector<std::tuple<uint, uint, uint, uint> >& coords;
    std::mutex& mtx;
    const std::vector<std::string>& seqs;
    bool print_status;

    void operator() (long i) const;

  HashSeq(param_struct& params_, const std::vector<std::string>& seqs_, btree::btree_map<uint, uint>& hashed_counts_, const std::vector<std::tuple<uint, uint,uint,uint> >& coords_, std::mutex& mtx_, bool print_status_) : params(params_), seqs(seqs_), hashed_counts(hashed_counts_), coords(coords_), mtx(mtx_), print_status(print_status_) {}
    

    ~HashSeq() {}
  };

};

