#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "SpookyV2.h"
#include <dlib/logger.h> 
#include <dlib/compress_stream.h>
#include <dlib/serialize.h>
#include "klib/kseq.h"
#include <zlib.h>
#include "cpp-btree/btree_map.h"
#include <mutex>
#include "data_set.h"

KSEQ_INIT(gzFile, gzread)

#define MAX_UINT_VAL 4294967295 //2^32-1
#define DP(i,j)     dp[(i) * (y.length() + 1) + (j)]
#define DPS(i,j)    dps[(i) * y.length() + (j)]

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
  std::string indices_filename;
  std::string ids_filename;
  bool is_ref;
  bool paired_end;
};


struct string_kernel2 {
  typedef float scalar_type;
  typedef std::string sample_type;

  typedef typename dlib::default_memory_manager mem_manager_type;

  string_kernel2(const int length, const double lambda,
		const string norm)
  : length(length), lambda(lambda), norm(norm)
  {}

  ~string_kernel2() {}

  const int length;
  const double lambda;
  const string norm;

  scalar_type operator() (const sample_type& x, const sample_type& y) const {
    float *dps, *dp;
    float kern[length];
    int i, j, l;


    /* Case a: both sequences empty */
    if (x.length() == 0 && y.length() == 0)
      return 1.0;

    /* Case b: one sequence empty */
    if (x.length() == 0 || y.length() == 0)
      return 0.0;

    /* Allocate temporary memory */
    dp = (float*) malloc(sizeof(float) * (x.length() + 1) * (y.length() + 1));
    dps = (float*) malloc(sizeof(float) * x.length() * y.length());
    if (!dp || !dps) {
      perror("Could not allocate memory for subsequence kernel");
      return 0;
    }

    /* Initalize dps */
    for (i = 0; i < x.length(); i++)
      for (j = 0; j < y.length(); j++) {
	//if (!strcmp(tolower(x.substr(i,1).c_str()),tolower(y.substr(j,1).c_str())))
	if (x[i] != y[j])
	  DPS(i,j) = lambda * lambda;
	else
	  DPS(i,j) = 0;
      }

    /* Initialize dp */
    for (i = 0; i < x.length() + 1; i++)
      DP(i, 0) = 0;
    for (j = 0; j < y.length() + 1; j++)
      DP(0, j) = 0;

    for (l = 0; l < length; l++) {
      kern[l] = 0;
      for (i = 0; i < x.length(); i++) {
	for (j = 0; j < y.length(); j++) {
	  DP(i + 1, j + 1) = DPS(i,j) + lambda * DP(i, j + 1) +
	    lambda * DP(i + 1, j) - lambda * lambda * DP(i, j);
	  //if (!strcmp(tolower(x.substr(i,1).c_str()),tolower(y.substr(j,1).c_str()))) {
	  if (x[i] != y[j]) {
	    kern[l] = kern[l] + DPS(i,j);
	    DPS(i,j) = lambda * lambda * DP(i, j);
	  }
	}
      }
    }
    
    free(dps);
    free(dp);
    
    return kern[length - 1];
  }

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
  bool count_hashed_kmers(param_struct&, btree::btree_map<std::string, btree::btree_map<uint, uint> >&, bool);
  bool hash_seq(std::string&, uint, btree::btree_map<uint, uint>&, std::mutex&);
  bool hash_seq(const std::vector<std::string>&, uint, btree::btree_map<std::string, btree::btree_map<uint, uint> >&, const std::string&, std::mutex&);
  template<typename T> bool add_dataset(const std::vector<T>&, const std::string&);
  bool add_taxonomy(const std::string&);
  uint hash_kmer(const std::string&);


  class BuilderTask {
    public:
      param_struct params;

      void execute();
      void check_memory();
      void hash_genome();
      void hash_sequences();
    
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

  class HashSeqs {
  public:
    param_struct& params;
    btree::btree_map<std::string, btree::btree_map<uint, uint> >& hashed_counts;
    const std::vector<std::tuple<std::vector<std::string>, std::string, uint> >& jobs;
    std::mutex& mtx;
    bool print_status;

    void operator() (long i) const;

  HashSeqs(param_struct& params_, btree::btree_map<std::string, btree::btree_map<uint, uint> >& hashed_counts_, const std::vector<std::tuple<std::vector<std::string>, std::string, uint> >& jobs_, std::mutex& mtx_, bool print_status_) : params(params_), hashed_counts(hashed_counts_), jobs(jobs_), mtx(mtx_), print_status(print_status_) {}


    ~HashSeqs() {}
  };


};

struct string_kernel {
  typedef double scalar_type;
  typedef std::string sample_type;

  typedef typename dlib::default_memory_manager mem_manager_type;

  string_kernel(const float c, const int normalize, const int symbol_size,
		const size_t max_length, int kn, double lambda)
  : _c(c), _normalize(normalize), _symbol_size(symbol_size),
    _max_length(max_length), _kn(kn), _lambda(lambda),
    _string_data(0), _kernel(0)
  {}

  string_kernel()
  : _c(0), _normalize(0), _symbol_size(0),
    _max_length(0), _kn(0), _lambda(0),
    _string_data(0), _kernel(0)
  {}

  ~string_kernel() {
    size_t size = (_string_data == 0) ? _size : _string_data->size();

    for (size_t i = 0; i < size; i++)
      delete[] _kernel[i];
    delete [] _kernel;
    delete _string_data;
  }

  /** Set the dataset to be used by the kernel. */
  void set_data(const std::vector<std::string> &strings) {
    assert(strings.size() > 0);
    _string_data = new DataSet(_max_length, _symbol_size);
    _string_data->load_strings(strings);
  }

  /** Calculate the kernel. */
  void compute_kernel() {
    assert(_string_data);

    // Initialize kernel 
    _kernel = new scalar_type *[_string_data->size()];
    for (size_t i = 0; i < _string_data->size(); i++)
      _kernel[i] = new scalar_type[_string_data->size()];

    // Start with all K filled with -1, then only calculate kernels as needed 
    for (size_t i = 0; i < _string_data->size(); i++)
      for (size_t j = 0; j < _string_data->size(); j++)
	_kernel[i][j] = -1;


    // Get values for normalization, it is computed for elements in diagonal 
    std::vector<scalar_type> norms(_string_data->size());
    if (_normalize) {
      for (size_t i = 0; i < _string_data->size(); i++) {
	norms[i] = kernel(_string_data->elements()[i], _string_data->elements()[i]);
	_kernel[i][i] = 1;
      }
    }

    // Compute kernel using dynamic programming 
    run_kernel_dp(norms, _kernel);
  }

  /** Return pointer to kernel matrix. */
  scalar_type **values() const {
    assert(_kernel);
    return _kernel;
  }

  void run_kernel_dp(const std::vector<scalar_type> &norms, scalar_type **K) const {
    assert(_string_data);
    for (size_t i = 0; i < _string_data->size(); i++) {
      for (size_t j = 0; j < _string_data->size(); j++) {
	if (K[j][i] == -1) {
	  K[j][i] = kernel(_string_data->elements()[j], _string_data->elements()[i]);
	  if (_normalize)
	    K[j][i] /= sqrt(norms[i] * norms[j]);
	  K[i][j] = K[j][i];
	}
      }
    }
  }


  /** Return the size of the NxN kernel. */
  size_t size() const {
    //assert(_string_data);
    return (_string_data == 0) ? _size : _string_data->size();
  }



  float _c;
  int _normalize;
  int _symbol_size;
  size_t _max_length;
  int _kn;
  double _lambda;
  DataSet *_string_data;
  scalar_type **_kernel;
  size_t _size;

  scalar_type operator() (const sample_type& a, const sample_type& b) const {
    // convert a -> DataElement x, b -> DataElement y
    DataElement x, y;
    size_t str_len = 0;
    std::string a_rep, b_rep;

    a_rep = (a < KMerge::rev_comp(a)) ? a : KMerge::rev_comp(a);
    b_rep = (b < KMerge::rev_comp(b)) ? b : KMerge::rev_comp(b);


    str_len = a_rep.length();
    x.allocate(str_len);
    for (size_t j = 0; j < str_len; j++) {
      char temp = tolower(*(a_rep.substr(j, 1).c_str()));
      assert(static_cast<int>(temp) < _symbol_size);
      x.attributes[j] = static_cast<int>(temp);
    }

    str_len = b_rep.length();
    y.allocate(str_len);
    for (size_t j = 0; j < str_len; j++) {
      char temp = tolower(*(b_rep.substr(j, 1).c_str()));
      assert(static_cast<int>(temp) < _symbol_size);
      y.attributes[j] = static_cast<int>(temp);
    }
    
    return kernel(x, y);
  }

  scalar_type kernel(const DataElement &x, const DataElement &y) const {
    scalar_type **Kd[2];
    for (int i = 0; i < 2; i++) {
      Kd[i] = new scalar_type *[x.length + 1];
      for (int j = 0; j < x.length + 1; j++) {
	Kd[i][j] = new scalar_type[y.length + 1];
      }
    }

    // Dynamic programming 
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < (x.length + 1); j++) {
	for (int k = 0; k < (y.length + 1); k++) {
	  Kd[i][j][k] = (i + 1) % 2;
	}
      }
    }

    // Calculate Kd and Kdd 
    for (int i = 1; i <= (_kn - 1); i++) {
      /* Set the Kd to zero for those lengths of s and t
	 where s (or t) has exactly length i-1 and t (or s)
	 has length >= i-1. L-shaped upside down matrix */
      for (int j = (i - 1); j <= (x.length - 1); j++) {
	Kd[i % 2][j][i - 1] = 0;
      }
      for (int j = (i - 1); j <= (y.length - 1); j++) {
	Kd[i % 2][i - 1][j] = 0;
      }
      for (int j = i; j <= (x.length - 1); j++) {
	scalar_type Kdd = 0;
	for (int m = i; m <= (y.length - 1); m++) {
	  if (x.attributes[j - 1] != y.attributes[m - 1]) {
	    // ((.))-1 is because indices start with 0 (not with 1)
	    Kdd = _lambda * Kdd;
	  } else {
	    Kdd = _lambda * (Kdd + (_lambda * Kd[(i + 1) % 2][j - 1][m - 1]));
	  }
	  Kd[i % 2][j][m] = _lambda * Kd[i % 2][j - 1][m] + Kdd;
	}
      }
    }

    // Calculate K 
    scalar_type sum = 0;
    for (int i = _kn; i <= x.length; i++) {
      for (int j = _kn; j <= y.length; j++) {
	if (x.attributes[((i)) - 1] == y.attributes[((j)) - 1]) {
	  sum += _lambda * _lambda * Kd[(_kn - 1) % 2][i - 1][j - 1];
	}
      }
    }
    
    // Delete 
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < x.length + 1; i++) {
	delete[] Kd[j][i];
      }
    }
    for (int i = 0; i < 2; i++) {
      delete[] Kd[i];
    }

    return sum;
  }

  bool operator== (
        const string_kernel& k
		   ) const 
  {
    if (_c != k._c)
      return false;

    if (_normalize != k._normalize)
      return false;

    if (_symbol_size != k._symbol_size)
      return false;

    if (_max_length != k._max_length)
      return false;

    if (_kn != k._kn)
      return false;

    if (_lambda != k._lambda)
      return false;

    if (_size != k._size)
      return false;

    for (size_t i = 0; i < _size; i++) {
      for (size_t j = 0; j < _size; j++) {
	if (_kernel[i][j] != k._kernel[i][j])
	  return false;
      }
    }

    // everything is equal
    return true;

  }


};

inline void serialize ( const string_kernel& item, std::ostream& out) {
  size_t size = (item._string_data == 0) ? item._size : item._string_data->size();
  dlib::serialize(size, out);
  for (size_t i = 0; i < size; i++) {
    for (size_t j = 0; j < size; j++) {
      dlib::serialize(item._kernel[i][j], out);
    }
  }
  
  // save the state of the kernel to the output stream
  dlib::serialize(item._c, out);
  dlib::serialize(item._normalize, out);
  dlib::serialize(item._symbol_size, out);
  dlib::serialize(item._max_length, out);
  dlib::serialize(item._kn, out);
  dlib::serialize(item._lambda, out);
}

inline void deserialize ( string_kernel& item, std::istream& in) {
  dlib::deserialize(item._size, in);
  item._kernel = new double*[item._size];
  for (size_t i = 0; i < item._size; i++) {
    item._kernel[i] = new double [item._size];
    for (size_t j = 0; j < item._size; j++) {
      dlib::deserialize(item._kernel[i][j], in);
    }
  }

  dlib::deserialize(item._c, in);
  dlib::deserialize(item._normalize, in);
  dlib::deserialize(item._symbol_size, in);
  dlib::deserialize(item._max_length, in);
  dlib::deserialize(item._kn, in);
  dlib::deserialize(item._lambda, in);
}
