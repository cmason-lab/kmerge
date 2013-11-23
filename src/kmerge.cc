#include "kmerge.h"
#include "lookup3.c"

using namespace std;

#define ROW_SIZE 0
#define COL_SIZE 1

KMerge::KMerge (const H5std_string& file_name) {
  this->file = new H5File( file_name, H5F_ACC_TRUNC );
}

uint KMerge::hashKmer(const std::string& kMer) {
  return((uint) hashlittle(kMer.c_str(), kMer.length(), 0));
}

bool KMerge::addDatasetToHDF5File(const H5std_string& ds_name, const hsize_t* chunk_dims, const hsize_t data_size, const uint* data) {
  /*                                                                                                                                                                                                                                        
   * Create property list for a dataset and set up fill values.                                                                                                                                                                             
   */

  uint fillvalue = 0;   /* Fill value for the dataset */
  DSetCreatPropList plist;
  hsize_t max_dims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};

  plist.setChunk(2, chunk_dims);
  plist.setFillValue(PredType::NATIVE_UINT, &fillvalue);
  plist.setDeflate( KMerge::GZIP_BEST_COMPRESSION );

  /*                                                                                                                                                                                                                                        
   * Create dataspace for the dataset in the file.                                                                                                                                                                                          
   */
  const hsize_t data_dims[2] = {data_size, 1}; 
  uint rank = (data_dims[ROW_SIZE] > data_dims[COL_SIZE]) ? data_dims[ROW_SIZE]:data_dims[COL_SIZE];
                                                                                                     
  DataSpace file_space( rank, data_dims, max_dims);

  /*                                                                                                                                                                                                                                        
   * Create dataset and write it into the file.                                                                                                                                                                                             
   */

  DataSpace mem_space( rank, data_dims );

  DataSet* dataset = new DataSet(file->createDataSet(ds_name, PredType::NATIVE_UINT, file_space, plist));
  dataset->write( data, PredType::NATIVE_UINT, mem_space, file_space );
  
  /*                                                                                                                                                                                                                                        
   * Close the dataset and free memory.
   */
  dataset->close();
  delete dataset;

  return true;
}

KMerge::~KMerge () {
  this->file->close();
  delete this->file;
}
