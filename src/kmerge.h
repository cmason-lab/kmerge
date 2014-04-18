#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <pthread.h>
#include "H5Cpp.h"
#include "fq.h"
#include "hdf5file.h"

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
using std::cout;
using std::endl;
using std::string;
#endif  // H5_NO_STD                                                                                                                                                                                                                          
#endif

using namespace std;

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

typedef unsigned int uint;

class KMerge;

struct param_struct {
  KMerge * kmerge;
  H5std_string hdf5_file_name;
  uint k_val_start;
  uint k_val_end;
  std::string seq_filename;
  std::string group_name;
  H5std_string hash_dataset_name;
  H5std_string counts_dataset_name;
} ;



class KMerge {
 private:
  H5std_string file_name;
  H5File *file;
  HDF5 *hdf5_file;
  static pthread_mutex_t mutex;

 public:
  static const uint GZIP_BEST_COMPRESSION = 9;
  static const uint CHUNK_ROW_SIZE = 1000000; //An estimate of the max number of kmer hashes produced for an organism


  KMerge(const H5std_string&);
  ~KMerge();
  //static uint hashKmer(const std::string&);
  static int hashKmer(const std::string&);
  static uint hash_kmer(const std::string&);
  static bool count_hashed_kmers(std::string&, uint, std::map<uint, uint>&);
  static bool add_hash_and_count(std::vector<uint>&, std::vector<uint>&, uint, uint);
  static bool add_hash_and_count(std::map<uint, uint>&, uint, uint);
  bool add_dataset(const std::string, uint, const uint*);
  static bool sort_kmer_hashes_and_counts(std::vector<uint>&, std::vector<uint>&);
  bool addHashAndCount(std::vector<uint>&, std::vector<uint>&, uint, uint);
  bool addHashAndCount(std::map<uint, uint>&, uint, uint);
  template <class T>
  std::vector<T> getDatasetFromHDF5File(const H5std_string&, const DataType&);
  bool addDatasetToHDF5File(const H5std_string&, const H5std_string&, const hsize_t, const uint*, const bool);
  static void addDatasetToHDF5FileT(void*);
  bool parseKmerCountsFile(const std::string&, std::vector<uint>&, std::vector<uint>&);
  bool parseKmerCountsFile(const std::string&, std::map<uint, uint>&);
  static void parseKmerCountsFileT(void*);
  static void parseAndWriteInThread(const param_struct&);
  bool appendToDataset(void);
  bool updateDataset(void);
  uint * getDatasetValues(void);
  std::vector<std::string> getGroupsFromHDF5File(const int&);

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

template <class T>
std::vector<T> KMerge::getDatasetFromHDF5File(const H5std_string& ds_name, const DataType& dt) {
  this->file = new H5File( this->file_name, H5F_ACC_RDONLY );

  /*                                                                                                                                                                                                                                                                                                                                                                       
   * Open the specified file and the specified dataset in the file.                                                                                                                                                                                                                                                                                                        
   */
  DataSet dataset = this->file->openDataSet( ds_name );

  /*                                                                                                                                                                                                                                                                                                                                                                       
   * Get filespace of the dataset.                                                                                                                                                                                                                                                                                                                                         
   */
  DataSpace file_space = dataset.getSpace();


  /*                                                                                                                                                                                                                                                                                                                                                                       
   * Get the number of dimensions in the dataspace.                                                                                                                                                                                                                                                                                                                        
   */
  int rank = file_space.getSimpleExtentNdims();


  /*                                                                                                                                                                                                                                                                                                                                                                       
   * Get the dimension size of each dimension in the dataspace                                                                                                                                                                                                                                                                                                             
   */
  hsize_t dims[2];
  int ndims = file_space.getSimpleExtentDims( dims, NULL);


  /*                                                                                                                                                                                                                                                                                                                                                                       
   * Define the memory space to read dataset.                                                                                                                                                                                                                                                                                                                              
   */

  DataSpace mem_space(rank, dims);


  /*                                                                                                                                                                                                                                                                                                                                                                    
   * Read dataset.                                                                                                                                                                                                                                                                                                                                                         
   */
  std::vector<T> data(dims[0]);
  dataset.read( &data[0], dt, mem_space, file_space );


  this->file->close();
  delete this->file;


  return data;
}
