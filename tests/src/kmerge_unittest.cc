#define CATCH_CONFIG_MAIN

#include <limits.h>
#include "H5Cpp.h"
#include "kmerge.h"
#include <math.h>
#include "queryProcessor.h"
#include "indexBuilder.h"
#include "catch.hpp"
#include "armadillo"

using namespace arma;

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
  
  KMerge* kmerge;

  HashTestFixture() {


    const string kmer1("ACTGA");
    const int kmer1_count = 4;
    const string kmer2("ATCGT");
    const int kmer2_count = 3;

    uint kmer1_hash_val = KMerge::hashKmer(kmer1);
    uint kmer2_hash_val = KMerge::hashKmer(kmer2);

    const H5std_string HDF5_FILE_NAME( "/home/darryl/Development/kmerge/tests/example.h5" );
    const H5std_string GROUP_NAME( "/org1" );
    const H5std_string HASH_DATASET_NAME( "/org1/kmer_hash" );
    const H5std_string COUNT_DATASET_NAME( "/org1/count" );

    try {
      this->kmerge = new KMerge(HDF5_FILE_NAME);
      std::vector<uint> hashes;
      std::vector<uint> counts;
                                                                                                                                 
      kmerge->addHashAndCount(hashes, counts, KMerge::hashKmer(kmer1), kmer1_count);
      kmerge->addHashAndCount(hashes, counts, KMerge::hashKmer(kmer2), kmer2_count);
      kmerge->addDatasetToHDF5File(GROUP_NAME, HASH_DATASET_NAME, hashes.size(), &hashes[0], true);
      kmerge->addDatasetToHDF5File(GROUP_NAME, COUNT_DATASET_NAME, counts.size(), &counts[0], false);

    } 
    catch( FileIException error ) {
      error.printError();
    }

    // catch failure caused by the DataSet operations
    catch( DataSetIException error ) {
	error.printError();
    }

    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error ) {
	error.printError();
    }

  }

  virtual ~HashTestFixture() {
    // Code here will be called immediately after each test (right
    // before the destructor).

    delete this->kmerge;

    if (remove( "/home/darryl/Development/kmerge/tests/example.h5" ) != 0) {
      perror( "Error deleting file");
    }
  }

  // Objects declared here can be used by all tests in the test case for Foo.
};

//TEST_CASE_METHOD(HashTestFixture, "HashKmerAndAppendToHDF5", "[HashTest]") {
  /* TODO:
   * New Implementation:
   * Extract hashes and counts for the group for which data is to be appended. Store each in a 
   * vector. Search the hashes vector for each hash to add. If its not there, add it at the end
   * and add count at the end of the counts vector. If it is there, add the count to the hashes
   * previous count. When all counts have been processed, write the data back to the HDF5 file
   * (either overwrite previous data or write new data and remove previous data).
   */


  //const string kmer1("ACTGA");
  //const int kmer1_count = 4;
  //const string kmer2("ATCGT");
  //const int kmer2_count = 3;
  //const string kmer3("ACTCA");
  //const int kmer3_count = 1;
  //const string kmer4("ATTTT");
  //const int kmer4_count = 2;

  //const H5std_string GROUP_NAME( "/org1" );
  //const H5std_string HASH_DATASET_NAME( "/org1/kmer_hash" );
  //const H5std_string COUNT_DATASET_NAME( "/org1/count" );
 

  //std::vector<uint> hashes = kmerge->getDatasetFromHDF5File(HASH_DATASET_NAME);
  //std::vector<uint> counts = kmerge->getDatasetFromHDF5File(COUNT_DATASET_NAME);

  //kmerge->addHashAndCount(hashes, counts, KMerge::hashKmer(kmer1), kmer1_count);
  //kmerge->addHashAndCount(hashes, counts, KMerge::hashKmer(kmer2), kmer2_count);
  //kmerge->addHashAndCount(hashes, counts, KMerge::hashKmer(kmer3), kmer3_count);
  //kmerge->addHashAndCount(hashes, counts, KMerge::hashKmer(kmer4), kmer4_count);
  /* Going to need to append data which has a hash already in the dataset and
   * add new hash for those which do not for any updates to the k-mer composition of a genome.
   * This is not neccessary when writing a dataset to the file the first time.
   */ 
  //kmerge->addDatasetToHDF5File(GROUP_NAME, HASH_DATASET_NAME, hashes.size(), &hashes[0], false);
  //kmerge->addDatasetToHDF5File(GROUP_NAME, COUNT_DATASET_NAME, counts.size(), &counts[0], false);


  //hashes.clear();
  //counts.clear();

  //REQUIRE(hashes.empty() == true);
  //REQUIRE(counts.empty() == true);

  //hashes = kmerge->getDatasetFromHDF5File(HASH_DATASET_NAME);
  //counts = kmerge->getDatasetFromHDF5File(COUNT_DATASET_NAME);

  //REQUIRE(hashes.size() == counts.size());
  //REQUIRE(hashes.size() == 4);

  //uint pos;
  //std::vector<uint>::iterator iter;

  //pos = 0;
  //for (std::vector<uint>::iterator iter=hashes.begin(); iter!=hashes.end(); iter++) {
  //  if (*iter==KMerge::hashKmer(kmer1)) {
  //    break;
  //  }
  //  pos++;
  //}

  //REQUIRE(iter != hashes.end());
  //REQUIRE(counts[pos] == kmer1_count);

  //pos = 0;
  //for (iter=hashes.begin(); iter!=hashes.end(); iter++) {
  //  if (*iter==KMerge::hashKmer(kmer2)) {
  //    break;
  //  }
  //  pos++;
  //}

  //REQUIRE(iter != hashes.end());
  //REQUIRE(counts[pos] == kmer2_count);

  //pos = 0;
  //for (iter=hashes.begin(); iter!=hashes.end(); iter++) {
  //  if (*iter==KMerge::hashKmer(kmer3)) {
  //    break;
  //  }
  //  pos++;
  //}

  //REQUIRE(iter != hashes.end());
  //REQUIRE(counts[pos] == kmer3_count);

  //pos = 0;
  //for (iter=hashes.begin(); iter!=hashes.end(); iter++) {
  //  if (*iter==KMerge::hashKmer(kmer4)) {
  //    break;
  //  }
  //  pos++;
  //}

  //REQUIRE(iter != hashes.end());
  //REQUIRE(counts[pos] == kmer4_count);

//}

TEST_CASE("ParseKmerCountsAndCreateHDF5", "[HashTest]") {
  std::vector<uint> hashes;
  std::vector<uint> counts;

  const std::string KMER_COUNT_FILE_NAME("/home/darryl/Development/kmerge/tests/sample_kmer_counts.txt");

  const H5std_string HDF5_FILE_NAME( "/home/darryl/Development/kmerge/tests/parse_example.h5" );
  const H5std_string GROUP_NAME( "/org1" );
  const H5std_string HASH_DATASET_NAME( "/org1/kmer_hash" );
  const H5std_string COUNT_DATASET_NAME( "/org1/count" );
  std::vector<uint>::iterator iter;


  try {
    KMerge* kmerge = new KMerge(HDF5_FILE_NAME);
  
    bool success = kmerge->parseKmerCountsFile(KMER_COUNT_FILE_NAME, hashes, counts);
    REQUIRE(success == true);

    REQUIRE(hashes.size() == 2);
    REQUIRE(counts.size() == 2);

    const string kmer1("ACTGA");
    const uint kmer1_count = 4;
    const string kmer2("ATCGT");
    const uint kmer2_count = 5;


    REQUIRE(hashes[0] == KMerge::hashKmer(kmer1));
    REQUIRE(hashes[1] == KMerge::hashKmer(kmer2));
    REQUIRE(counts[0] == kmer1_count);
    REQUIRE(counts[1] == kmer2_count);
  
    kmerge->addDatasetToHDF5File(GROUP_NAME, HASH_DATASET_NAME, hashes.size(), &hashes[0], true);
    kmerge->addDatasetToHDF5File(GROUP_NAME, COUNT_DATASET_NAME, counts.size(), &counts[0], false);


    std::vector<uint> hashes2 = kmerge->getDatasetFromHDF5File(HASH_DATASET_NAME);                                                                                                                        
    std::vector<uint> counts2 = kmerge->getDatasetFromHDF5File(COUNT_DATASET_NAME);

    
  
    REQUIRE(hashes2.size() == 2);
    REQUIRE(counts2.size() == 2);


    uint pos = 0; 
    for (iter=hashes2.begin(); iter!=hashes2.end(); iter++) {
      if (*iter==KMerge::hashKmer(kmer1)) {      
	break;
      }                                                                                                                                                                                               
      pos++;                                                                                                                                                              
    }
    REQUIRE(iter != hashes2.end());
    REQUIRE(counts2[pos] == kmer1_count);

    pos = 0;
    
    for (iter=hashes2.begin(); iter!=hashes2.end(); iter++) {
      if (*iter==KMerge::hashKmer(kmer2)) {
	break;
      }
      pos++;
    }

    REQUIRE(iter != hashes2.end());
    REQUIRE(counts2[pos] == kmer2_count);


    if (remove("/home/darryl/Development/kmerge/tests/parse_example.h5" ) != 0) {
      perror( "Error deleting file");
    }

    
  }    
  catch( FileIException error ) {
    error.printError();
    cout << "File error" << std::endl;
    if (remove("/home/darryl/Development/kmerge/tests/parse_example.h5" ) != 0) {
      perror( "Error deleting file");
    }
  }

  // catch failure caused by the DataSet operations 
  catch( DataSetIException error ) {
    error.printError();
    cout << "Dataset error" << std::endl;
    if (remove("/home/darryl/Development/kmerge/tests/parse_example.h5" ) != 0) {
      perror( "Error deleting file");
    }
  }

  // catch failure caused by the DataSpace operations 
  catch( DataSpaceIException error ) {
    error.printError();
    cout << "Dataspace error" << std::endl;
    if (remove("/home/darryl/Development/kmerge/tests/parse_example.h5" ) != 0) {
      perror( "Error deleting file");
    }
  }
}


TEST_CASE("WriteSparseMatrixToHDF5AndReadFromAramdillo", "[SparseHDF5Test]") {
  const string kmer1("ACTGA");
  const uint kmer1_count = 4;
  const string kmer2("ATCGT");
  const uint kmer2_count = 3;

  uint kmer1_hash_val = KMerge::hashKmer(kmer1);
  uint kmer2_hash_val = KMerge::hashKmer(kmer2);

  const H5std_string FILE_NAME( "/home/darryl/Development/kmerge/tests/sparse.h5" );
  const H5std_string GROUP_NAME( "/org1" );
  const H5std_string HASH_DATASET_NAME( "/org1/kmer_hash" );
  const H5std_string COUNT_DATASET_NAME( "/org1/count" );
  hsize_t maxdims = H5S_UNLIMITED;
  hsize_t chunk_dims = 2;
  const int FSPACE_RANK = 1;
  const int FSPACE_DIM = 2;

  /*
   * Create a file and group;
   */

  H5File* file = new H5File( FILE_NAME, H5F_ACC_TRUNC );

  Group* group = new Group( file->createGroup( "/org1" ));

  /*                                                                                                                             
   * Create property list for a dataset and set up fill values.                                                                           
   */

  uint fillvalue = 0;   /* Fill value for the dataset */
  DSetCreatPropList plist;
  plist.setChunk(1, &chunk_dims);
  plist.setFillValue(PredType::NATIVE_UINT, &fillvalue);
  plist.setDeflate( KMerge::GZIP_BEST_COMPRESSION );

  /*                                                                                                 
   * Create dataspace for the dataset in the file.  
   */

                                       
  hsize_t fdim = FSPACE_DIM; // dim size of ds (on disk)
  DataSpace fspace( FSPACE_RANK, &fdim, &maxdims);

  /*                               
   * Create datasets and write it into the file.           
   */


  DataSet* hash_dataset = new DataSet(file->createDataSet(HASH_DATASET_NAME, PredType::NATIVE_UINT, fspace, plist));
  DataSet* count_dataset = new DataSet(file->createDataSet(COUNT_DATASET_NAME, PredType::NATIVE_UINT, fspace, plist));

  /* 
   * Create dataspace for the datasets.
   */

  hsize_t dim1 = FSPACE_DIM;

  /* 
   * Dimension size of the first dataset (in memory)  
   */

  DataSpace mspace1( FSPACE_RANK, &dim1 );
  uint vector[FSPACE_DIM];     

  // vector buffer for dset                                                                                                                                                                        
  vector[0] = kmer1_hash_val;
  vector[1] = kmer2_hash_val;

  hash_dataset->write( vector, PredType::NATIVE_UINT, mspace1, fspace );

  vector[0] = kmer1_count;
  vector[1] = kmer2_count;

  count_dataset->write( vector, PredType::NATIVE_UINT, mspace1, fspace );

  /*  
   * Close the dataset and the file.     
   */


  hash_dataset->close();
  delete hash_dataset;
  count_dataset->close();
  delete count_dataset;
  group->close();
  delete group;
  file->close();
  delete file;

  /*                                                                                                                                       
   * Open the specified file and the specified dataset in the file.                                                                            
   */
  H5File in_file( FILE_NAME, H5F_ACC_RDONLY );
  DataSet in_hash_dataset = in_file.openDataSet( HASH_DATASET_NAME );
  DataSet in_count_dataset = in_file.openDataSet( COUNT_DATASET_NAME );
  /*                                                                                                    
   * Get dataspace of the dataset.
   */
  DataSpace filespace = in_hash_dataset.getSpace();

  /*                                                                                                                   
   * Get the number of dimensions in the dataspace. 
   */
  int hash_rank = filespace.getSimpleExtentNdims();

  REQUIRE(hash_rank == 1);

  /*                                                                                          
   * Get the dimension size of each dimension in the dataspace and                                                                        
   * display them.                                                                                                            
   */
  hsize_t dims;
  int ndims = filespace.getSimpleExtentDims( &dims, NULL);
  
  REQUIRE(dims == (hsize_t) 2);
  REQUIRE(ndims == 1);

  /*                                                                                                          
   * Define the memory space to read dataset.                                                                                                       
   */
  DataSpace mspace(hash_rank, &dims);

  /*                                                                                                                                                
   * Read dataset.                                                                                                                 
   */
  std::vector<uword> hashed_kmers(dims);
  in_hash_dataset.read( &hashed_kmers[0], PredType::NATIVE_ULLONG, mspace, filespace );


  REQUIRE(kmer1_hash_val == hashed_kmers[0]);
  REQUIRE(kmer2_hash_val == hashed_kmers[1]);

  uint *hash_counts = new uint[dims];
  in_count_dataset.read( hash_counts, PredType::NATIVE_UINT, mspace, filespace );

  REQUIRE(kmer1_count == hash_counts[0]);
  REQUIRE(kmer2_count == hash_counts[1]);


  uword temp[] = {0, dims}; // need to have 1 more entry than the number of columns which is greater than the length of "counts"
  uvec org_start(temp, dims);

  uvec hashes(hashed_kmers);
  Col<uint> counts(hash_counts, dims);
  SpMat<uint> A(hashes, org_start, counts, UINT_MAX, org_start.n_elem - 1);

  REQUIRE(A(kmer1_hash_val, 0) == kmer1_count);
  REQUIRE(A(kmer2_hash_val, 0) == kmer2_count);

  delete [] hash_counts;
}
