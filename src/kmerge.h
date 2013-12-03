#include <string>
#include "H5Cpp.h"
#include "fq.h"

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
using std::cout;
using std::endl;
using std::string;
#endif  // H5_NO_STD                                                                                                                                                                                                                          
#endif

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

using namespace std;
typedef unsigned int uint;

class KMerge {
 private:
  H5std_string file_name;
  H5File *file;

 public:
  static const uint GZIP_BEST_COMPRESSION = 9;
  static const uint CHUNK_ROW_SIZE = 1000000; //An estimate of the max number of kmer hashes produced for an organism
  
  KMerge(const H5std_string&);
  ~KMerge();
  static uint hashKmer(const std::string&);
  bool addHashAndCount(std::vector<uint>&, std::vector<uint>&, uint, uint);
  std::vector<uint> getDatasetFromHDF5File(const H5std_string&);
  bool addDatasetToHDF5File(const H5std_string&, const H5std_string&, const hsize_t, const uint*, const bool);
  bool parseKmerCountsFile(const std::string&, std::vector<uint>&, std::vector<uint>&);
  bool appendToDataset(void);
  bool updateDataset(void);
  uint * getDatasetValues(void);
};
