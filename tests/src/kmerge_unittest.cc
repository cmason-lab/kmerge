#define CATCH_CONFIG_MAIN

#include <limits.h>
#include "H5Cpp.h"
#include "kmerge.h"
#include <math.h>
#include "queryProcessor.h"
#include "indexBuilder.h"
#include "catch.hpp"
#include "armadillo"
#include "threadpool.h"

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

TEST_CASE_METHOD(HashTestFixture, "AddAnotherOrganism", "[HashTest]") {
  std::vector<uint> hashes;
  std::vector<uint> counts;

  const std::string KMER_COUNT_FILE_NAME("/home/darryl/Development/kmerge/tests/k3.counts.gz");

  const H5std_string GROUP_NAME( "/org2" );
  const H5std_string HASH_DATASET_NAME( "/org2/kmer_hash" );
  const H5std_string COUNT_DATASET_NAME( "/org2/count" );
  std::vector<uint>::iterator iter;

  try {

    bool success = this->kmerge->parseKmerCountsFile(KMER_COUNT_FILE_NAME, hashes, counts);
    REQUIRE(success == true);

    REQUIRE(hashes.size() == 64);
    REQUIRE(counts.size() == 64);

    const string kmer0("AAA");
    const string kmer63("TTT");
    const uint kmer0_count = 95944;
    const uint kmer63_count = 99230;

    REQUIRE(hashes[0] == KMerge::hashKmer(kmer0));
    REQUIRE(hashes[63] == KMerge::hashKmer(kmer63));
    REQUIRE(counts[0] == kmer0_count);
    REQUIRE(counts[63] == kmer63_count);

    this->kmerge->addDatasetToHDF5File(GROUP_NAME, HASH_DATASET_NAME, hashes.size(), &hashes[0], true);
    this->kmerge->addDatasetToHDF5File(GROUP_NAME, COUNT_DATASET_NAME, counts.size(), &counts[0], false);


    std::vector<uint> hashes2 = this->kmerge->getDatasetFromHDF5File(HASH_DATASET_NAME);
    std::vector<uint> counts2 = this->kmerge->getDatasetFromHDF5File(COUNT_DATASET_NAME);



    REQUIRE(hashes2.size() == 64);
    REQUIRE(counts2.size() == 64);


    uint pos = 0;
    for (iter=hashes2.begin(); iter!=hashes2.end(); iter++) {
      if (*iter==KMerge::hashKmer(kmer0)) {
        break;
      }
      pos++;
    }
    REQUIRE(iter != hashes2.end());
    REQUIRE(counts2[pos] == kmer0_count);

    pos = 0;

    for (iter=hashes2.begin(); iter!=hashes2.end(); iter++) {
      if (*iter==KMerge::hashKmer(kmer63)) {
        break;
      }
      pos++;
    }

    REQUIRE(iter != hashes2.end());
    REQUIRE(counts2[pos] == kmer63_count);

    REQUIRE(counts == counts2);
    REQUIRE(hashes == hashes2);


  }
  catch( FileIException error ) {
    error.printError();
    cout << "File error" << std::endl;
  }

  // catch failure caused by the DataSet operations                                                                                                                                                                  
  catch( DataSetIException error ) {
    error.printError();
    cout << "Dataset error" << std::endl;
  }

  // catch failure caused by the DataSpace operations                                                                                                                                                                
  catch( DataSpaceIException error ) {
    error.printError();
    cout << "Dataspace error" << std::endl;
  }


}

TEST_CASE("ParseKmerCountsAndCreateHDF5", "[HashTest]") {
  std::vector<uint> hashes;
  std::vector<uint> counts;
  std::map<uint, uint> hashed_counts;

  const std::string KMER_COUNT_FILE_NAME("/home/darryl/Development/kmerge/tests/k3.counts.gz");

  const H5std_string HDF5_FILE_NAME( "/home/darryl/Development/kmerge/tests/parse_example.h5" );
  const H5std_string GROUP_NAME( "/org1" );
  const H5std_string HASH_DATASET_NAME( "/org1/kmer_hash" );
  const H5std_string COUNT_DATASET_NAME( "/org1/count" );
  std::vector<uint>::iterator iter;


  try {
    KMerge* kmerge = new KMerge(HDF5_FILE_NAME);
  
    bool success = kmerge->parseKmerCountsFile(KMER_COUNT_FILE_NAME, hashed_counts);
    REQUIRE(success == true);

    const string kmer0("AAA");
    const string kmer63("TTT");
    uint kmer0_pos = 0;
    uint kmer63_pos = 0;
    uint pos = 0;

    for (std::map<uint, uint>::iterator map_iter = hashed_counts.begin(); map_iter != hashed_counts.end(); map_iter++) {
      if (map_iter->first == KMerge::hashKmer(kmer0)) {
	kmer0_pos = pos;
      }
      if (map_iter->first == KMerge::hashKmer(kmer63)) {
	kmer63_pos = pos;
      }
      hashes.push_back(map_iter->first);
      counts.push_back(map_iter->second);
      pos++;
    }

    REQUIRE(hashes.size() == 64);
    REQUIRE(counts.size() == 64);


    const uint kmer0_count = 95944;
    const uint kmer63_count = 99230;


    REQUIRE(hashes[kmer0_pos] == KMerge::hashKmer(kmer0));
    REQUIRE(hashes[kmer63_pos] == KMerge::hashKmer(kmer63));
    REQUIRE(counts[kmer0_pos] == kmer0_count);
    REQUIRE(counts[kmer63_pos] == kmer63_count);
    
    kmerge->addDatasetToHDF5File(GROUP_NAME, HASH_DATASET_NAME, hashes.size(), &hashes[0], true);
    kmerge->addDatasetToHDF5File(GROUP_NAME, COUNT_DATASET_NAME, counts.size(), &counts[0], false);


    std::vector<uint> hashes2 = kmerge->getDatasetFromHDF5File(HASH_DATASET_NAME);                                                                                                                        
    std::vector<uint> counts2 = kmerge->getDatasetFromHDF5File(COUNT_DATASET_NAME);

    
  
    REQUIRE(hashes2.size() == 64);
    REQUIRE(counts2.size() == 64);

    
    pos = 0; 
    for (iter=hashes2.begin(); iter!=hashes2.end(); iter++) {
      if (*iter==KMerge::hashKmer(kmer0)) {      
	break;
      }                                                                                                                                                                                               
      pos++;                                                                                                                                                              
    }
    REQUIRE(iter != hashes2.end());
    REQUIRE(counts2[pos] == kmer0_count);

    pos = 0;
    
    for (iter=hashes2.begin(); iter!=hashes2.end(); iter++) {
      if (*iter==KMerge::hashKmer(kmer63)) {
	break;
      }
      pos++;
    }

    REQUIRE(iter != hashes2.end());
    REQUIRE(counts2[pos] == kmer63_count);

    REQUIRE(counts == counts2);
    REQUIRE(hashes == hashes2);

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

TEST_CASE("ThreadedParseKmerCountsAndCreateHDF5", "[HashTest]") {
  std::vector<uint> hashes;
  std::vector<uint> counts;

  param_struct params1;
  param_struct params2;
  param_struct params3;

  params1.k_val_start = 3;
  params1.k_val_end = 3;
  params2.k_val_start = 3;
  params2.k_val_end = 3;
  params3.k_val_start = 3;
  params3.k_val_end = 3;


  params1.hdf5_file_name =  "/home/darryl/Development/kmerge/tests/thread_example.h5";
  params2.hdf5_file_name = params1.hdf5_file_name;
  params3.hdf5_file_name = params2.hdf5_file_name;

  params1.group_name =  "15660";
  params2.group_name = "165199";
  params3.group_name ="29309";

  params1.hash_dataset_name =  "/15660/kmer_hash";
  params1.count_dataset_name = "/15660/count";

  params2.hash_dataset_name = "/165199/kmer_hash";
  params2.count_dataset_name = "/165199/count";

  params3.hash_dataset_name = "/29309/kmer_hash";
  params3.count_dataset_name = "/29309/count";

  uint thread_count = 3;

  ThreadPool tp(thread_count);


  try {
    KMerge* kmerge = new KMerge(params1.hdf5_file_name);

    params1.kmerge = kmerge;
    params2.kmerge = kmerge;
    params3.kmerge = kmerge;


    Task t1(&KMerge::parseAndWriteInThread, (void*) &params1);
    Task t2(&KMerge::parseAndWriteInThread, (void*) &params2);
    Task t3(&KMerge::parseAndWriteInThread, (void*) &params3);

    tp.initialize_threadpool();
    

    tp.add_task(&t1);
    tp.add_task(&t2);
    tp.add_task(&t3);

    sleep(2);

    tp.destroy_threadpool();


    std::vector<uint> hashes1 = kmerge->getDatasetFromHDF5File("/15660/kmer_hash");
    std::vector<uint> hashes2 = kmerge->getDatasetFromHDF5File("/165199/kmer_hash");
    std::vector<uint> hashes3 = kmerge->getDatasetFromHDF5File("/29309/kmer_hash");

    std::vector<uint> counts1 = kmerge->getDatasetFromHDF5File("/15660/count");
    std::vector<uint> counts2 = kmerge->getDatasetFromHDF5File("/165199/count");
    std::vector<uint> counts3 = kmerge->getDatasetFromHDF5File("/29309/count");

    const string kmer0("AAA");
    const string kmer63("TTT");
    uint kmer0_count = 85060;
    uint kmer63_count = 85174;
    uint kmer0_pos = 0;
    uint kmer63_pos = 0;
    uint pos = 0;

    for (std::vector<uint>::iterator iter = hashes1.begin(); iter != hashes1.end(); iter++) {
      if ((*iter) == KMerge::hashKmer(kmer0)) {
	kmer0_pos = pos;
      }
      if ((*iter) == KMerge::hashKmer(kmer63)) {
	kmer63_pos = pos;
      }
      pos++;
    }

    REQUIRE(hashes1[kmer0_pos] == KMerge::hashKmer(kmer0));
    REQUIRE(hashes1[kmer63_pos] == KMerge::hashKmer(kmer63));
    REQUIRE(counts1[kmer0_pos] == kmer0_count);
    REQUIRE(counts1[kmer63_pos] == kmer63_count);


    kmer0_count = 65798;
    kmer63_count = 65803;
    pos = 0;
    for (std::vector<uint>::iterator iter = hashes2.begin(); iter != hashes2.end(); iter++) {
      if ((*iter) == KMerge::hashKmer(kmer0)) {
        kmer0_pos = pos;
      }
      if ((*iter) == KMerge::hashKmer(kmer63)) {
	kmer63_pos = pos;
      }
      pos++;
    }

    REQUIRE(hashes2[kmer0_pos] == KMerge::hashKmer(kmer0));
    REQUIRE(hashes2[kmer63_pos] == KMerge::hashKmer(kmer63));
    REQUIRE(counts2[kmer0_pos] == kmer0_count);
    REQUIRE(counts2[kmer63_pos] == kmer63_count);


    kmer0_count = 95944;
    kmer63_count = 99230;

    pos = 0;
    for (std::vector<uint>::iterator iter = hashes3.begin(); iter != hashes3.end(); iter++) {
      if ((*iter) == KMerge::hashKmer(kmer0)) {
        kmer0_pos = pos;
      }
      if ((*iter) == KMerge::hashKmer(kmer63)) {
	kmer63_pos = pos;
      }
      pos++;
    }

    REQUIRE(hashes3[kmer0_pos] == KMerge::hashKmer(kmer0));
    REQUIRE(hashes3[kmer63_pos] == KMerge::hashKmer(kmer63));
    REQUIRE(counts3[kmer0_pos] == kmer0_count);
    REQUIRE(counts3[kmer63_pos] == kmer63_count);


    if (remove("/home/darryl/Development/kmerge/tests/thread_example.h5" ) != 0) {
      perror( "Error deleting file");
    }

    delete kmerge;
  }
  catch( FileIException error ) {
    error.printError();
    cout << "File error" << std::endl;
    if (remove("/home/darryl/Development/kmerge/tests/thread_example.h5" ) != 0) {
      perror( "Error deleting file");
    }
  }
  
  // catch failure caused by the DataSet operations                                                                          
  catch( DataSetIException error ) {
    error.printError();
    cout << "Dataset error" << std::endl;
    if (remove("/home/darryl/Development/kmerge/tests/thread_example.h5" ) != 0) {
      perror( "Error deleting file");
    }
  }
  
  // catch failure caused by the DataSpace operations
  catch( DataSpaceIException error ) {
    error.printError();
    cout << "Dataspace error" << std::endl;
    if (remove("/home/darryl/Development/kmerge/tests/thread_example.h5" ) != 0) {
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
