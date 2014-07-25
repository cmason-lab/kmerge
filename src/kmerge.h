#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <pthread.h>
#include "fq.h"
#include "hdf5file.h"
#include "SpookyV2.h"
#include <dlib/logger.h> 
#include "klib/kseq.h"
#include <zlib.h>
#include "cpp-btree/btree_map.h"

KSEQ_INIT(gzFile, gzread)

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
  uint num_threads;
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
  KMerge(const std::string&, const std::string&, const std::string&);
  ~KMerge();
  static std::string rev_comp(const std::string&);
  void build(param_struct&);
  bool count_hashed_kmers_fasta(param_struct&, btree::btree_map<uint, uint>&);
  bool count_hashed_kmers_fastq(param_struct&, btree::btree_map<uint, uint>&);
  bool add_dataset(const std::string, uint, const uint*, HDF5*);
  bool add_taxonomy(const std::string&);
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

  class CountAndHashSeqFasta {
  public:

    param_struct& params;
    btree::btree_map<uint, uint>& hashed_counts;
    std::string seq;
    pthread_mutex_t& m;

    void operator() (long i) const;

  CountAndHashSeqFasta( param_struct& params_, btree::btree_map<uint, uint>& hashed_counts_, const std::string seq_, pthread_mutex_t& m_) : params(params_), hashed_counts(hashed_counts_), seq(seq_), m(m_) {}

    ~CountAndHashSeqFasta() {}
    };

  class CountAndHashSeqFastq {
  public:
    
    param_struct& params;
    btree::btree_map<uint, uint>& hashed_counts;
    pthread_mutex_t& m;
    
    void operator() (long i) const;

  CountAndHashSeqFastq( param_struct& params_, btree::btree_map<uint, uint>& hashed_counts_, pthread_mutex_t& m_ ) : params(params_), hashed_counts(hashed_counts_), m(m_) {}

    ~CountAndHashSeqFastq() {}
  };

};

