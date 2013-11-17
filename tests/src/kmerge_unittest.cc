#define CATCH_CONFIG_MAIN

#include <limits.h>
#include "H5Cpp.h"
#include "kmerge.h"
#include <math.h>
#include "queryProcessor.h"
#include "indexBuilder.h"
#include "catch.hpp"


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


// The fixture for testing class k-mer hashing.
class HashTestFixture  {
 public:
  // You can remove any or all of the following functions if its body
  // is empty.

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:
  

  HashTestFixture() {
    // Code here will be called immediately after the constructor (right
    // before each test).

    const string kmer1("ACTGA");
    const int kmer1_count = 4;
    const string kmer2("ATCGT");
    const int kmer2_count = 3;

    uint kmer1_hash_val = hashKmer(kmer1);
    uint kmer2_hash_val = hashKmer(kmer2);

    const H5std_string FILE_NAME( "/home/darryl/Development/kmerge/tests/example.h5" );
    const H5std_string DATASET_NAME( "test_ds" );
    hsize_t maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
    hsize_t chunk_dims[2] = {1, 2};
    const int FSPACE_RANK = 2;
    const int FSPACE_DIM1 = 2;
    const int FSPACE_DIM2 = 2;

    /*
     * Create a file.
     */

    H5File* file = new H5File( FILE_NAME, H5F_ACC_TRUNC );

    /*                                                                                                                             
     * Create property list for a dataset and set up fill values.                                                                           
     */

    uint fillvalue = 0;   /* Fill value for the dataset */
    DSetCreatPropList plist;
    plist.setChunk(2, chunk_dims);
    plist.setFillValue(PredType::NATIVE_UINT, &fillvalue);

    /*                                                                                                 
     * Create dataspace for the dataset in the file.  
     */

    hsize_t fdim[] = {FSPACE_DIM1, FSPACE_DIM2}; // dim sizes of ds (on disk)                                       
    DataSpace fspace( FSPACE_RANK, fdim, maxdims);

    /*                               
     * Create dataset and write it into the file.           
     */


    DataSet* dataset = new DataSet(file->createDataSet(DATASET_NAME, PredType::NATIVE_UINT, fspace, plist));

    /* 
     * Create dataspace for the first dataset.
     */

    hsize_t dim1[] = {FSPACE_DIM1, FSPACE_DIM2};

    /* 
     * Dimension size of the first dataset (in memory)  
     */

    DataSpace mspace1( FSPACE_RANK, dim1 );
    uint vector[FSPACE_DIM1][FSPACE_DIM2];     

    // vector buffer for dset                                                                                                                                                                        
    vector[0][0] = kmer1_hash_val;
    vector[0][1] = kmer1_count;
    vector[1][0] = kmer2_hash_val;
    vector[1][1] = kmer2_count;

    dataset->write( vector, PredType::NATIVE_UINT, mspace1, fspace );

    /*  
     * Close the dataset and the file.     
     */

    delete dataset;
    delete file;

  }

  virtual ~HashTestFixture() {
    // Code here will be called immediately after each test (right
    // before the destructor).

    //if (remove( "/home/darryl/Development/kmerge/tests/example.h5" ) != 0) {
    //  perror( "Error deleting file");
    //}
  }

  // Objects declared here can be used by all tests in the test case for Foo.
};


TEST_CASE_METHOD(HashTestFixture, "HashKmerAndWriteToHDF5", "[HashTest]") {
  const string kmer1("ACTGA");
  const int kmer1_count = 4;
  const string kmer2("ATCGT");
  const int kmer2_count = 3;

  uint kmer1_hash_val = hashKmer(kmer1);
  uint kmer2_hash_val = hashKmer(kmer2);
  REQUIRE(612248635 == kmer1_hash_val);
  REQUIRE(2960106173 == kmer2_hash_val);
 
  const H5std_string FILE_NAME( "/home/darryl/Development/kmerge/tests/example.h5" );
  const H5std_string DATASET_NAME( "test_ds" );
  //hsize_t maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
  //hsize_t chunk_dims[2] = {1, 2};
  //const int FSPACE_RANK = 2;
  //const int FSPACE_DIM1 = 2;
  //const int FSPACE_DIM2 = 2;
  

  /*                                                                                                                                       
   * Open the specified file and the specified dataset in the file.                                                                            
   */
  H5File in_file( FILE_NAME, H5F_ACC_RDONLY );
  DataSet in_dataset = in_file.openDataSet( DATASET_NAME );

  /*                                                                                                    
   * Get dataspace of the dataset.
   */
  DataSpace in_filespace = in_dataset.getSpace();

  /*                                                                                                                   
   * Get the number of dimensions in the dataspace. 
   */
  int rank = in_filespace.getSimpleExtentNdims();

  /*                                                                                          
   * Get the dimension size of each dimension in the dataspace and                                                                        
   * display them.                                                                                                            
   */
  hsize_t dims[2];
  int ndims = in_filespace.getSimpleExtentDims( dims, NULL);
  
  /*                                                                                                          
   * Define the memory space to read dataset.                                                                                                       
   */
  DataSpace mspace2(rank, dims);

  /*                                                                                                                                                
   * Read dataset back and display.                                                                                                                 
   */
  int data_out[dims[0]][dims[1]];  // buffer for dataset to be read                                        
  in_dataset.read( data_out, PredType::NATIVE_UINT, mspace2, in_filespace );

  REQUIRE(kmer1_hash_val == data_out[0][0]);
  REQUIRE(kmer1_count == data_out[0][1]);
  REQUIRE(kmer2_hash_val == data_out[1][0]);
  REQUIRE(kmer2_count == data_out[1][1]);
}

TEST_CASE_METHOD(HashTestFixture, "HashKmerAndAppendToHDF5", "[HashTest]") {
  /*      
   * Open the specified file and the specified dataset in the file.
   */
  const H5std_string FILE_NAME( "/home/darryl/Development/kmerge/tests/example.h5" );
  const H5std_string DATASET_NAME( "test_ds" );

  H5File in_file( FILE_NAME, H5F_ACC_RDWR );
  DataSet in_dataset = in_file.openDataSet( DATASET_NAME );
  
  int nlines = 1;
  int ncols = 2;
  hsize_t dimsext[2] = {nlines, ncols};
  uint vector[nlines][ncols];
  const string kmer_to_append("ACTGA");
  uint kmer_to_append_hash_val = hashKmer(kmer_to_append);
  vector[0][0] = kmer_to_append_hash_val;
  uint kmer_to_append_count = 2;
  vector[0][1] = kmer_to_append_count;
  const string kmer1("ACTGA");
  const int kmer1_count = 4;
  const string kmer2("ATCGT");
  const int kmer2_count = 3;

  uint kmer1_hash_val = hashKmer(kmer1);
  uint kmer2_hash_val = hashKmer(kmer2);

  /*
   * Get dimensions from the dataset to append to and create memory space.
   */
  hsize_t dims1[2];
  DataSpace file_space = in_dataset.getSpace();
  int ndims = file_space.getSimpleExtentDims( dims1, NULL);
  int rank = file_space.getSimpleExtentNdims();

  /*
   * Extend dataset to allow for new data.
   */
  hsize_t size[2];
  size[0] = dims1[0] + dimsext[0];
  size[1] = dims1[1];
  in_dataset.extend(size);

  /*
   * Select hyperslab on file dataset.
   */
  DataSpace *file_space2 = new DataSpace( in_dataset.getSpace() );
  hsize_t offset[2] = {dims1[0], 0};
  file_space2->selectHyperslab(H5S_SELECT_SET, dimsext, offset);

  /*
   * Define memory space.
   */

  DataSpace *mem_space = new DataSpace(2, dimsext, NULL);


  /*
   * Write to the extended portion of the dataset.
   */
  in_dataset.write( vector, PredType::NATIVE_UINT, *mem_space, *file_space2);

  delete file_space2;
  delete mem_space;
  in_file.close();
  in_dataset.close();

  /*                                                   
   * Read dataset back and display.                                                                  
   */

  const H5std_string FILE_NAME2( "/home/darryl/Development/kmerge/tests/example.h5" );
  const H5std_string DATASET_NAME2( "test_ds" );

  H5File in_file2( FILE_NAME2, H5F_ACC_RDONLY );
  DataSet in_dataset2 = in_file2.openDataSet( DATASET_NAME2 );

  DataSpace file_space3 = in_dataset2.getSpace();
  DataSpace mem_space2(rank, size);


  int data_out[size[0]][size[1]];  // buffer for dataset to be read                                                                                                                                                                          
  in_dataset2.read( data_out, PredType::NATIVE_UINT, mem_space2, file_space3 );

  REQUIRE(kmer1_hash_val == data_out[0][0]);
  REQUIRE(kmer1_count == data_out[0][1]);
  REQUIRE(kmer2_hash_val == data_out[1][0]);
  REQUIRE(kmer2_count == data_out[1][1]);
  REQUIRE(kmer_to_append_hash_val == data_out[2][0]);
  REQUIRE(kmer_to_append_count == data_out[2][1]);

}

TEST_CASE_METHOD(HashTestFixture, "HashKmerAndFindVallue", "[HashTest]") {
  const H5std_string FILE_NAME( "/home/darryl/Development/kmerge/tests/example.h5" );
  const H5std_string DATASET_NAME( "test_ds" );
  const string kmer("ATCGT");
  uint hashed_kmer = hashKmer(kmer);
  const std::vector<double> hash_val (1, hashed_kmer);
  std::vector<uint64_t> coords;

  /* 
   * Build indexes on example file.
   */
  IndexBuilder *index_builder = new IndexBuilder(FILE_NAME, FQ::FQ_HDF5, "", 3);
  index_builder->buildIndexes();

  delete index_builder;

  /*
   * Initialize query processor and search for kmer1 value.
   */

  QueryProcessor* query_processor = new QueryProcessor(FILE_NAME, FQ::FQ_HDF5);
  int num_found = query_processor->executeEqualitySelectionQuery(DATASET_NAME, hash_val, coords);
  

  /*
   * Ensure that value is found and get the selected value.
   */
  REQUIRE(num_found == 1);
  uint64_t* data = new uint64_t[1];
  bool found = query_processor->getSelectedData(DATASET_NAME, coords, data);
  

  /*
   * Make sure value is what we expect it to be.
   */
  REQUIRE(found == 1);
  REQUIRE(data[0] == hashed_kmer);
  coords[1]++;
  found = query_processor->getSelectedData(DATASET_NAME, coords, data);
  REQUIRE(found == 1);
  REQUIRE(data[0] == 3);
  
  delete data;
  delete query_processor;


  /*
   * Update value from 3 and change to 5.
   */

  H5File in_file( FILE_NAME, H5F_ACC_RDWR );
  DataSet in_dataset = in_file.openDataSet( DATASET_NAME );

  
  /*
   * Select hyperslab on file dataset.
   */
  hsize_t single_row[2] = {1, 2};
  DataSpace *file_space = new DataSpace( in_dataset.getSpace() );
  hsize_t offset[2] = {coords[0], 0};
  file_space->selectHyperslab(H5S_SELECT_SET, single_row, offset);

  // Overwrite data at found coordinates.
  DataSpace *mem_space = new DataSpace(2, single_row, NULL);
  uint vector[1][2];
  vector[0][0] = hashed_kmer;
  uint updated_kmer_count = 5;
  vector[0][1] = updated_kmer_count;
  in_dataset.write( vector, PredType::NATIVE_UINT, *mem_space, *file_space);

  in_file.close();

  /*
   * Make sure kmer count at that position is now 5.
   */

  QueryProcessor* query_processor2 = new QueryProcessor(FILE_NAME, FQ::FQ_HDF5);

  uint64_t data2;
  found = query_processor2->getSelectedData(DATASET_NAME, coords, &data2);
  REQUIRE(found == 1);
  REQUIRE(data2 == 5);
}

