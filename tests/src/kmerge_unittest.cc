#define CATCH_CONFIG_MAIN

#include <limits.h>
#include "kmerge.h"
#include <math.h>
#include "hdf5file.h"
#include "catch.hpp"
#include <sstream>
#include <dlib/threads.h>
#include <dlib/serialize.h>
#include <dlib/logger.h>
#include "cpp-btree/btree_map.h"
#include <iomanip>
#include "indexBuilder.h"
#include "queryProcessor.h"
#include "H5Cpp.h"
#include "blosc_filter.h"


#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

TEST_CASE("WriteMatrixToHDF5File", "HDF5Test") {
  uint num_cols = 1;
  unsigned int cd_values[7];
  char *version, *date;
  int r, i;
  const H5std_string ds_name("counts");
  const H5std_string file_name("matrix.h5");

  /* Register the filter with the library */
  r = register_blosc(&version, &date);
  printf("Blosc version info: %s (%s)\n", version, date);


  H5::H5File *file = new H5::H5File( file_name, H5F_ACC_TRUNC );
  
/*ifstream ifile("matrix.h5");
  if (ifile) {
    this->file = new H5File( this->file_name, H5F_ACC_RDWR );
  } else {
    this->file = new H5File( this->file_name, H5F_ACC_TRUNC );
    }*/

  /*                              
   * Create property list for a dataset and set up fill values. 
   */
  uint total_columns = 3125;
  uint rank = 2;
  uint fillvalue = 0;   /* Fill value for the dataset */
  H5::DSetCreatPropList plist;
  hsize_t max_dims[2] = {MAX_UINT_VAL, H5S_UNLIMITED};
  hsize_t chunk_dims[2] = {KMerge::CHUNK_ROW_SIZE, 1};

  cd_values[4] = 1;       /* compression level */
  cd_values[5] = 1;       /* 0: shuffle not active, 1: shuffle active */
  cd_values[6] = BLOSC_LZ4HC; /* the actual compressor to use */

  plist.setChunk(rank, chunk_dims);
  plist.setFillValue(H5::PredType::NATIVE_UINT, &fillvalue);
  plist.setFilter(FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values);
  plist.setShuffle();
  plist.setFletcher32();

  /* 
   * Create dataspace for the dataset                   
   */
  const hsize_t data_dims[2] = {MAX_UINT_VAL, num_cols}; 
                                                                                                     
  H5::DataSpace file_space( rank, data_dims, max_dims);



  /* 
   * Create dataset and write it into the file.                        
   */

  H5::DataSpace mem_space( rank, data_dims, max_dims );

  H5::DataSet* dataset = new H5::DataSet( file->createDataSet(ds_name, H5::PredType::NATIVE_UINT, file_space, plist));

  uint hash;
  std::stringstream path;
  btree::btree_map<uint, uint> hashed_counts;
  param_struct params;

  params.k_val_start = 5;
  params.k_val_end = 5;
  params.hdf5_filename = "dummy.h5";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/genome.test.fasta.gz";
  params.tmp_hashes_filename = "";
  params.tmp_counts_filename = "";
  params.group_name = "/test";
  params.hash_dataset_name =  "/test/kmer_hash";
  params.counts_dataset_name = "/test/count";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  // build sample kmerge file                                                                                                                                                                                                                      

  KMerge *kmerge = new KMerge(params.hdf5_filename.c_str(), "lookup3", ".");
  params.kmerge = kmerge;
  bool success = kmerge->count_hashed_kmers_fasta(params, hashed_counts);
  REQUIRE(success == true);

  REQUIRE(hashed_counts.size() == 324);

  std::vector<uint> counts(MAX_UINT_VAL);

  std::cout << "loading counts from btree_map" << std::endl;
  for (btree::btree_map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    counts[iter->first] = iter->second;
  }
  std::cout << "Done" << std::endl;

  hashed_counts.clear();
  btree::btree_map<uint, uint>().swap( hashed_counts );

  std::cout << "Writing data" << std::endl;
  dataset->write( &counts[0], H5::PredType::NATIVE_UINT, mem_space, file_space );
  std::cout << "Done" << std::endl;

  counts.clear();
  std::vector<uint>().swap(counts);

  dataset->close();
  delete dataset;

  file->close();
  delete file;


  // search for data                                                                                                                                                                                                                               
  QueryProcessor* query_processor = new QueryProcessor("matrix.h5", FQ::FQ_HDF5);

  // test that search does not return incorrect values                                                                                                                                                                                             

  //find coordinates of values of interest                                                                                                                                                                                                         
  uint hashes_arr[] = {1022408118, 0, 1043881921, 0, 1045475391, 0, 1049435108, 0, 1130773752, 0};
  //uint hashes_arr[] = {1022408118, 1043881921, 1045475391, 1049435108, 1130773752};
  // test that search does not return incorrect values                                                                                                                                                                                             
  std::vector<uint64_t> coords(hashes_arr, hashes_arr + 10);
  //std::vector<uint64_t> coords(hashes_arr, hashes_arr + 5);
  uint counts_arr[] = {1, 1, 1, 4, 97};
  std::vector<uint> freq(counts_arr, counts_arr + 5);

  coords.push_back(1); // count = 0
  coords.push_back(0); // count = 0

  coords.push_back(10); // count = 0
  coords.push_back(0); // count = 0

  //use returned coordinates from above to get hashes from "kmer_hash" and counts from "count"
  uint num_hashes = coords.size() / 2;

  uint* c_result = new uint[num_hashes];
  query_processor->getSelectedData("counts", coords, c_result, "/");

  delete query_processor;

  uint adj_pos = 0;
  for(uint i = 0; i < num_hashes; i++) {
    if (c_result[i] != 0) {
      REQUIRE(freq[adj_pos] == c_result[i]);
      adj_pos++;
    }
  }

  delete [] c_result;

  file = new H5::H5File( file_name, H5F_ACC_RDWR );
  
  dataset = new H5::DataSet( file->openDataSet( ds_name ) );
  
  H5::DataSpace dataspace = dataset->getSpace();
  rank = dataspace.getSimpleExtentNdims();

  hsize_t dims_out[2];
  int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
  REQUIRE(dims_out[0] == MAX_UINT_VAL);
  REQUIRE(dims_out[1] == 1);
  hsize_t offset[2];
  offset[0] = 0;
  offset[1] = dims_out[1];
  dims_out[1]++; // extend dataset by 1
  dataset->extend( dims_out );

  DataSpace fspace = dataset->getSpace();


  fspace.selectHyperslab( H5S_SELECT_SET, data_dims, offset );
  DataSpace mspace( rank, data_dims );

  // build sample kmerge file

  kmerge = new KMerge(params.hdf5_filename.c_str(), "lookup3", ".");
  params.kmerge = kmerge;
  success = kmerge->count_hashed_kmers_fasta(params, hashed_counts);
  REQUIRE(success == true);

  REQUIRE(hashed_counts.size() == 324);

  counts.resize(MAX_UINT_VAL);

  for (btree::btree_map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    counts[iter->first] = iter->second;
  }

  hashed_counts.clear();
  btree::btree_map<uint, uint>().swap( hashed_counts );

  dataset->write( &counts[0], H5::PredType::NATIVE_UINT, mspace, fspace );

  counts.clear();
  std::vector<uint>().swap(counts);

  dataset->close();
  delete dataset;

  file->close();
  delete file;

  // search for data
  query_processor = new QueryProcessor("matrix.h5", FQ::FQ_HDF5);

  // test that search does not return incorrect values
  //find coordinates of values of interest 
  uint hashes_arr2[] = {1022408118, 0, 1043881921, 0, 1045475391, 0, 1049435108, 0, 1130773752, 0, 1022408118, 1, 1043881921, 1, 1045475391, 1, 1049435108, 1, 1130773752, 1};

  std::vector<uint64_t> coords2(hashes_arr2, hashes_arr2 + 20);

  //use returned coordinates from above to get hashes from "kmer_hash" and counts from "count"
  num_hashes = coords2.size() / 2;

  c_result = new uint[num_hashes];
  query_processor->getSelectedData("counts", coords2, c_result, "/");

  delete query_processor;


  adj_pos = 0;
  for(uint i = 0; i < num_hashes; i++) {
    if (c_result[i] != 0) {
      REQUIRE(freq[adj_pos % 5] == c_result[i]);
      adj_pos++;
    }
  }

  delete [] c_result;

  if (remove("matrix.h5") != 0) {
    perror( "Error deleting file");
  }
}

TEST_CASE("StoreKmerHashesAsFastBitIndices", "[FastBitTest]") {
  uint hash;
  std::stringstream path;
  btree::btree_map<uint, uint> hashed_counts;
  param_struct params;

  params.k_val_start = 5;
  params.k_val_end = 5;
  params.hdf5_filename = "index_data.h5";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/genome.test.fasta.gz";
  params.tmp_hashes_filename = "";
  params.tmp_counts_filename = "";
  params.group_name = "/test";
  params.hash_dataset_name =  "/test/kmer_hash";
  params.counts_dataset_name = "/test/count";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
 
  // build sample kmerge file

  KMerge *kmerge = new KMerge(params.hdf5_filename.c_str(), "lookup3", ".");
  params.kmerge = kmerge;
  bool success = kmerge->count_hashed_kmers_fasta(params, hashed_counts);
  REQUIRE(success == true);

  REQUIRE(hashed_counts.size() == 324);

  std::vector<uint> hashes, counts;
  for (btree::btree_map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    hashes.push_back(iter->first);                                                                
    counts.push_back(iter->second);
  }

  REQUIRE(params.kmerge->add_dataset(params.hash_dataset_name, hashes.size(), &hashes[0], NULL) == true);
  REQUIRE(params.kmerge->add_dataset(params.counts_dataset_name, counts.size(), &counts[0], NULL) == true);

  delete kmerge;
 
  // add indices to file
  IndexBuilder* index_builder = new IndexBuilder("index_data.h5", FQ::FQ_HDF5, "index_data.h5");
  index_builder->buildIndexes(0, "kmer_hash");
  delete index_builder;

  // search for data
  QueryProcessor* query_processor = new QueryProcessor("index_data.h5", FQ::FQ_HDF5);

  // test that search does not return incorrect values
  std::vector<double> np_ids;
  np_ids.push_back(10);
  std::vector<uint64_t> coords;
  uint hits = query_processor->executeEqualitySelectionQuery("kmer_hash", np_ids, coords, "/test");
  REQUIRE(hits == 0);
  REQUIRE(coords.size() == 0);

  //find coordinates of values of interest
  uint hashes_arr[] = {1022408118, 1043881921, 1045475391, 1049435108, 1130773752};
  std::vector<double> ids(hashes_arr, hashes_arr + 5);
  std::vector<uint> hashes_ref(hashes_arr, hashes_arr + 5);
  uint counts_arr[] = {1, 1, 1, 4, 97};
  std::vector<uint> freq(counts_arr, counts_arr + 5);

  ids.push_back(1); // count = 0
  ids.push_back(10); // count = 0

  hits = query_processor->executeEqualitySelectionQuery("kmer_hash", ids, coords, "/test");
  REQUIRE(hits == 5);
  REQUIRE(ids.size() != coords.size());
  REQUIRE(coords.size() == 5);

  //use returned coordinates from above to get hashes from "kmer_hash" and counts from "count"

  uint* h_result = new uint[coords.size()]; 
  uint* c_result = new uint[coords.size()];
  query_processor->getSelectedData("kmer_hash", coords, h_result, "/test");
  query_processor->getSelectedData("count", coords, c_result, "/test");

  delete query_processor;
  
  for(uint i = 0; i < coords.size(); i++) {
    std::vector<uint>::iterator iter = std::find (hashes_ref.begin(), hashes_ref.end(), h_result[i]);
    REQUIRE(iter != hashes_ref.end());
    uint pos = iter - hashes_ref.begin();
    REQUIRE(freq[pos] == c_result[i]);
  }
  
  delete [] h_result;
  delete [] c_result;


  if (remove("index_data.h5") != 0) {
    perror( "Error deleting file");
  }
 
}

TEST_CASE("StoreDenseCounts", "[HashTest]") {
  uint hash;
  std::stringstream path;
  btree::btree_map<uint, uint> hashed_counts;
  param_struct params;

  params.k_val_start = 5;
  params.k_val_end = 5;
  params.hdf5_filename = "dense_data.h5";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/genome.test.fasta.gz";
  params.tmp_hashes_filename = "";
  params.tmp_counts_filename = "";
  params.group_name = "/test";
  params.hash_dataset_name =  "/test/kmer_hash";
  params.counts_dataset_name = "/test/count";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
 
  // build sample kmerge file

  KMerge *kmerge = new KMerge(params.hdf5_filename.c_str(), "lookup3", ".");
  params.kmerge = kmerge;
  bool success = kmerge->count_hashed_kmers_fasta(params, hashed_counts);
  REQUIRE(success == true);

  REQUIRE(hashed_counts.size() == 324);

  std::vector<uint> counts(MAX_UINT_VAL);
  
  for (btree::btree_map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    counts[iter->first] = iter->second; 
  }

  REQUIRE(params.kmerge->add_dataset(params.counts_dataset_name, MAX_UINT_VAL, &counts[0], NULL) == true);

  delete kmerge;


  // search for data
  QueryProcessor* query_processor = new QueryProcessor("dense_data.h5", FQ::FQ_HDF5);

  // test that search does not return incorrect values

  //find coordinates of values of interest
  uint hashes_arr[] = {1022408118, 1043881921, 1045475391, 1049435108, 1130773752};
  // test that search does not return incorrect values
  std::vector<uint64_t> coords(hashes_arr, hashes_arr + 5);
  uint counts_arr[] = {1, 1, 1, 4, 97};
  std::vector<uint> freq(counts_arr, counts_arr + 5);

  coords.push_back(1); // count = 0
  coords.push_back(10); // count = 0

  //use returned coordinates from above to get hashes from "kmer_hash" and counts from "count"
 
  uint* c_result = new uint[coords.size()];
  query_processor->getSelectedData("count", coords, c_result, "/test");

  delete query_processor;
  uint adj_pos = 0;
  for(uint i = 0; i < coords.size(); i++) {
    if (c_result[i] != 0) {
      REQUIRE(freq[adj_pos] == c_result[i]);
      adj_pos++;
    }
  }
  
  delete [] c_result;


  if (remove("dense_data.h5") != 0) {
    perror( "Error deleting file");
  }
 
}

TEST_CASE("IndicesAreUsableFromMergedHDF5File", "[FastBitTest]") {
  uint hits;
  std::vector<uint64_t> coords;
  // search for data
  QueryProcessor* query_processor = new QueryProcessor("/home/darryl/Development/kmerge/tests/reference.u.h5", FQ::FQ_HDF5);


  //find coordinates of values of interest
  uint hashes_arr[] = {313855374, 883124211, 1275167949, 164111441, 291025028, 313855374, 2149088619, 1700750264, 2469142716, 2837084308};
  std::vector<double> ids(hashes_arr, hashes_arr + 10);
  std::vector<uint> hashes_ref(hashes_arr, hashes_arr + 10);
  uint counts1[] = {41872, 30947, 23680};
  std::vector<uint> freq1(counts1, counts1 + 3);

  hits = query_processor->executeEqualitySelectionQuery("kmer_hash", ids, coords, "/1");
  REQUIRE(hits == 3);
  REQUIRE(ids.size() != coords.size());
  REQUIRE(coords.size() == 3);


  //use returned coordinates from above to get hashes from "kmer_hash" and counts from "count"

  uint* h_result = new uint[coords.size()]; 
  uint* c_result = new uint[coords.size()];
  query_processor->getSelectedData("kmer_hash", coords, h_result, "/1");
  query_processor->getSelectedData("count", coords, c_result, "/1");
  
  for(uint i = 0; i < coords.size(); i++) {
    std::vector<uint>::iterator iter = std::find (hashes_ref.begin(), hashes_ref.end(), h_result[i]);
    REQUIRE(iter != hashes_ref.end());
    iter = std::find(freq1.begin(), freq1.end(), c_result[i]);
    REQUIRE(iter != freq1.end());
  }

  coords.clear();
  REQUIRE(coords.size() == 0);
  
  delete [] h_result;
  delete [] c_result;


  uint counts2[] = {44677, 36661, 41872, 30213, 22141};
  std::vector<uint> freq2(counts2, counts2 + 5);

  hits = query_processor->executeEqualitySelectionQuery("kmer_hash", ids, coords, "/2");
  REQUIRE(hits == 5);
  REQUIRE(ids.size() != coords.size());
  REQUIRE(coords.size() == 5);


  //use returned coordinates from above to get hashes from "kmer_hash" and counts from "count"

  h_result = new uint[coords.size()]; 
  c_result = new uint[coords.size()];
  query_processor->getSelectedData("kmer_hash", coords, h_result, "/2");
  query_processor->getSelectedData("count", coords, c_result, "/2");
  
  for(uint i = 0; i < coords.size(); i++) {
    std::vector<uint>::iterator iter = std::find (hashes_ref.begin(), hashes_ref.end(), h_result[i]);
    REQUIRE(iter != hashes_ref.end());
    iter = std::find(freq2.begin(), freq2.end(), c_result[i]);
    REQUIRE(iter != freq2.end());
  }

  coords.clear();
  REQUIRE(coords.size() == 0);
  
  delete [] h_result;
  delete [] c_result;


  uint counts3[] = {44677, 36661, 30143, 30213, 13127, 22141};
  std::vector<uint> freq3(counts3, counts3 + 6);

  hits = query_processor->executeEqualitySelectionQuery("kmer_hash", ids, coords, "/3");
  REQUIRE(hits == 6);
  REQUIRE(ids.size() != coords.size());
  REQUIRE(coords.size() == 6);


  //use returned coordinates from above to get hashes from "kmer_hash" and counts from "count"

  h_result = new uint[coords.size()]; 
  c_result = new uint[coords.size()];
  query_processor->getSelectedData("kmer_hash", coords, h_result, "/3");
  query_processor->getSelectedData("count", coords, c_result, "/3");
  
  for(uint i = 0; i < coords.size(); i++) {
    std::vector<uint>::iterator iter = std::find (hashes_ref.begin(), hashes_ref.end(), h_result[i]);
    REQUIRE(iter != hashes_ref.end());
    iter = std::find(freq3.begin(), freq3.end(), c_result[i]);
    REQUIRE(iter != freq3.end());
  }

  coords.clear();
  REQUIRE(coords.size() == 0);
  
  delete [] h_result;
  delete [] c_result;


  delete query_processor;
}

TEST_CASE("BTreeMapIsSortedTest", "[ContainerTest]") {
  btree::btree_map<uint, uint> m;
  uint prev = 0;

  m[1] = 1;
  m[4] = 4;
  m[7] = 7;

  REQUIRE(m.size() == 3);

  REQUIRE(m[1] == 1);
  REQUIRE(m[4] == 4);
  REQUIRE(m[7] == 7);

  for (btree::btree_map<uint, uint>::iterator m_iter = m.begin(); m_iter != m.end(); m_iter++) {
    REQUIRE(m_iter->first > prev);
    prev = m_iter->first;
  }

}

TEST_CASE("CountKmersTest", "[HashTest]") {

  int k = 5, l;
  std::map<std::string, uint> kmer_counts, jf_counts;
  kseq_t *seq;
  gzFile fp;

  fp = gzopen("/home/darryl/Development/kmerge/tests/genome.test.fasta.gz", "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s);
    for (int i = 0; i < seq->seq.l - k + 1; i++) {
      std::string kmer = seq_str.substr(i, k);
      std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
      if(kmer.find_first_not_of("ACGT") != std::string::npos) { // skip kmers containing non-nucleotides
	continue;
      }
      if (kmer_counts.find(kmer) != kmer_counts.end()) {
	kmer_counts[kmer]++;
      } else {
	kmer_counts[kmer] = 1;
      }
    }
  }
  kseq_destroy(seq);
  gzclose(fp);

  REQUIRE(kmer_counts.size() == 424);


  // compare to jellyfish results
  // perl ../scripts/contiguous_fasta.pl <(zcat genome.test.fa.gz) | ~/bin/fastx_toolkit/bin/fasta_formatter -w 80 | gzip > genome.test.contig.fa.gz
  // ~/bin/jellyfish count -m 5 -o output -s 10000 <(zcat genome.test.contig.fa.gz)
  // ~/bin/jellyfish dump -c output > genome.test.kmers.txt
  
  ifstream in_file ("/home/darryl/Development/kmerge/tests/genome.test.kmers.txt");

  std::string seq2;
  uint count;
  if (in_file.is_open()) {
    while ( !in_file.eof() ) {
      in_file >> seq2 >> count;
      if (seq2 != "") {
	jf_counts[seq2] = count;
	REQUIRE(kmer_counts[seq2] == count);
      }
    }
    in_file.close();
  }
  REQUIRE(jf_counts.size() == kmer_counts.size());
}

TEST_CASE("CountKmersInFastqFile", "[HashTest]") {

  int k = 7, l;
  std::map<std::string, uint> kmer_counts, jf_counts;
  kseq_t *seq;
  gzFile fp;

  fp = gzopen("/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz", "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s);
    for (int i = 0; i < seq->seq.l - k + 1; i++) {
      std::string kmer = seq_str.substr(i, k);
      std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
      if(kmer.find_first_not_of("ACGT") != std::string::npos) { // skip kmers containing non-nucleotides
	continue;
      }
      if (kmer_counts.find(kmer) != kmer_counts.end()) {
	kmer_counts[kmer]++;
      } else {
	kmer_counts[kmer] = 1;
      }
    }
  }
  kseq_destroy(seq);
  gzclose(fp);

  //~/bin/kanalyze-0.9.5/count -d 20 -f fastqgz -k 5 -o sample/sample.k5.txt sample/sample.fastq.gz

  ifstream in_file ("/home/darryl/Development/kmerge/tests/sample/sample.k7.txt");

  std::string seq2;
  uint count;
  if (in_file.is_open()) {
    while ( !in_file.eof() ) {
      in_file >> seq2 >> count;
      if (seq2 != "") {
	jf_counts[seq2] = count;
	REQUIRE(kmer_counts[seq2] == count);
      }
    }
    in_file.close();
  }
  REQUIRE(jf_counts.size() == kmer_counts.size());
}

TEST_CASE("CountHashedKmersInFastaFile", "[HashTest]") {
  btree::btree_map<uint, uint> hashed_counts, jf_hashed_counts;
  param_struct params;
  std::vector<std::string> files;

  params.k_val_start = 3;
  params.k_val_end = 7;
  params.hdf5_filename = "/home/darryl/Development/kmerge/tests/fasta.h5";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/208831/208831.fasta.gz";
  params.tmp_hashes_filename = "/home/darryl/Development/kmerge/tests/208831/hashes.bin";
  params.tmp_counts_filename = "/home/darryl/Development/kmerge/tests/208831/counts.bin";
  params.group_name = "/208831";
  params.hash_dataset_name =  "/208831/kmer_hash";
  params.counts_dataset_name = "/208831/count";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;  


  //~/bin/kanalyze-0.9.5/count -d 20 -f fastagz -k 3 -o 208831/sample.k3.txt 208831/208831.fasta.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastagz -k 5 -o 208831/sample.k5.txt 208831/208831.fasta.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastagz -k 7 -o 208831/sample.k7.txt 208831/208831.fasta.gz

  files.push_back("/home/darryl/Development/kmerge/tests/208831/sample.k3.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/208831/sample.k5.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/208831/sample.k7.txt");

  KMerge *kmerge = new KMerge(params.hdf5_filename.c_str(), "lookup3", ".");

  params.kmerge = kmerge;


  bool success = kmerge->count_hashed_kmers_fasta(params, hashed_counts);
  REQUIRE(success == true);
  

  std::string seq, last = "";
  uint count;
  for (std::string filename : files ) {
    ifstream in_file (filename.c_str());
    if (in_file.is_open()) {
      while ( !in_file.eof() ) {
	in_file >> seq >> count;
	if (seq != last) {
	  uint hash = kmerge->hash_kmer(seq);
	  jf_hashed_counts[hash] = jf_hashed_counts[hash] +  count;
	}
	last = seq;
      }
      in_file.close();
    }
  }

  delete kmerge;

  REQUIRE(jf_hashed_counts.size() == hashed_counts.size());
  
  for (btree::btree_map<uint, uint>::const_iterator b_it = jf_hashed_counts.begin(); b_it != jf_hashed_counts.end(); b_it++) {
    REQUIRE(hashed_counts[b_it->first] == b_it->second);
  }

  if (remove(params.hdf5_filename.c_str()) != 0) {
    perror( "Error deleting file");
  }
}

TEST_CASE("CountHashedKmersInFastqFile", "[HashTest]") {
  btree::btree_map<uint, uint> hashed_counts, jf_hashed_counts;
  param_struct params;
  std::vector<std::string> files;

  params.k_val_start = 3;
  params.k_val_end = 7;
  params.hdf5_filename = "fastq.h5";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz";
  params.tmp_hashes_filename = "/home/darryl/Development/kmerge/tests/sample/hashes.bin";
  params.tmp_counts_filename = "/home/darryl/Development/kmerge/tests/sample/counts.bin";
  params.group_name = "/sample";
  params.hash_dataset_name =  "/sample/kmer_hash";
  params.counts_dataset_name = "/sample/count";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;  

  //~/bin/kanalyze-0.9.5/count -d 20 -f fastqgz -k 3 -o sample/sample.k3.txt sample/sample.fastq.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastqgz -k 5 -o sample/sample.k5.txt sample/sample.fastq.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastqgz -k 7 -o sample/sample.k7.txt sample/sample.fastq.gz

  files.push_back("/home/darryl/Development/kmerge/tests/sample/sample.k3.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/sample/sample.k5.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/sample/sample.k7.txt");

  KMerge *kmerge = new KMerge(params.hdf5_filename.c_str(), "lookup3", ".");

  params.kmerge = kmerge;


  bool success = kmerge->count_hashed_kmers_fastq(params, hashed_counts);
  REQUIRE(success == true);
  

  std::string seq, last = "";
  uint count;
  for (std::string filename : files ) {
    ifstream in_file (filename.c_str());
    if (in_file.is_open()) {
      while ( !in_file.eof() ) {
	in_file >> seq >> count;
	if (seq != last) {
	  uint hash = kmerge->hash_kmer(seq);
	  jf_hashed_counts[hash] = jf_hashed_counts[hash] +  count;
	}
	last = seq;
      }
      in_file.close();
    }
  }

  delete kmerge;

  REQUIRE(jf_hashed_counts.size() == hashed_counts.size());

  for (btree::btree_map<uint, uint>::const_iterator b_it = jf_hashed_counts.begin(); b_it != jf_hashed_counts.end(); b_it++) {
    REQUIRE(hashed_counts[b_it->first] == b_it->second);
  }

  if (remove(params.hdf5_filename.c_str()) != 0) {
    perror( "Error deleting file");
  }
}

TEST_CASE("CountHashedKmersParallelFastq", "[HashTest]") {

  int l;
  btree::btree_map<uint, uint> hashed_counts, jf_hashed_counts;
  pthread_mutex_t mutex;
  KMerge *kmerge = new KMerge("dummy.h5", "lookup3", ".");
  param_struct params;
  std::vector<std::string> files;

  params.k_val_start = 3;
  params.k_val_end = 7;
  params.hdf5_filename = "";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz";
  params.tmp_hashes_filename = "";
  params.tmp_counts_filename = "";
  params.group_name = "/sample";
  params.hash_dataset_name =  "/sample/kmer_hash";
  params.counts_dataset_name = "/sample/count";
  params.kmerge = kmerge;
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;


  pthread_mutex_init(&mutex, NULL);
  KMerge::CountAndHashSeqFastq func(params, hashed_counts, mutex);
  dlib::parallel_for(params.num_threads, (params.k_val_start + 1) / 2, (params.k_val_end + 1) / 2 + 1, func);

  pthread_mutex_destroy(&mutex);
  REQUIRE(hashed_counts.size() == 8727);

  delete kmerge;
}


TEST_CASE("CountHashedKmersParallelFasta", "[HashTest]") {

  int l;
  btree::btree_map<uint, uint> hashed_counts;
  kseq_t *seq;
  gzFile fp;
  pthread_mutex_t mutex;
  KMerge *kmerge = new KMerge("dummy.h5", "lookup3", ".");
  param_struct params;

  params.k_val_start = 3;
  params.k_val_end = 11;
  params.hdf5_filename = "problem.h5";
  params.seq_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/208831/208831.fasta.gz";
  params.tmp_hashes_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/208831/hashes.bin";
  params.tmp_counts_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/208831/counts.bin";
  params.group_name = "/208831";
  params.hash_dataset_name =  "/208831/kmer_hash";
  params.counts_dataset_name = "/208831/count";
  params.kmerge = kmerge;
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;


  pthread_mutex_init(&mutex, NULL);
  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s);
    KMerge::CountAndHashSeqFasta func(params, hashed_counts, seq_str, mutex);
    dlib::parallel_for(params.num_threads, (params.k_val_start + 1) / 2, (params.k_val_end + 1) / 2 + 1, func);
  }

  delete kmerge;
  pthread_mutex_destroy(&mutex);
  REQUIRE(hashed_counts.size() == 1591922);
  kseq_destroy(seq);
  gzclose(fp);

}


TEST_CASE("TestHashedKmersAndReverseComplementReturnSameHashVal", "[HashTest]") {
  std::string test_seq("ACTAG");
  std::string rc_test_seq = KMerge::rev_comp(test_seq);
  
  REQUIRE(rc_test_seq == "CTAGT");
  KMerge *kmerge = new KMerge("dummy.h5", "lookup3", ".");
  
  uint hash1 = kmerge->hash_kmer(test_seq), hash2 = kmerge->hash_kmer(rc_test_seq);
  REQUIRE(hash1 == hash2);
  delete kmerge;
}

TEST_CASE("SerializeAndDeserializeVectors", "[SerializeTest]") {
  std::map<uint, uint> m;
  std::vector<uint> hashes, counts;
  
  m[0] = 4;
  m[1] = 1;
  m[2] = 0;
  m[3] = 2;

  for (std::map<uint, uint>::const_iterator m_iter = m.begin(); m_iter != m.end(); m_iter++) {
    hashes.push_back(m_iter->first);
    counts.push_back(m_iter->second);
  }

  ofstream out_hashes_file("tmp_hashes.bin", ios::out | ios::binary), out_counts_file("tmp_counts.bin", ios::out | ios::binary);
  try {
    dlib::serialize(hashes, out_hashes_file);
    dlib::serialize(counts, out_counts_file);
  } catch (dlib::serialization_error& e) {
    cout << "Unable to serialize data" << endl;
  }

  std::vector<uint>().swap( hashes );
  std::vector<uint>().swap( counts );

  out_hashes_file.close();
  out_counts_file.close();

  ifstream in_hashes_file("tmp_hashes.bin", ios::in | ios::binary), in_counts_file("tmp_counts.bin", ios::in | ios::binary);

  try {
    dlib::deserialize(hashes, in_hashes_file);
    dlib::deserialize(counts, in_counts_file);
  } catch (dlib::serialization_error& e) {
    cout << "Unable to deserialize data" << endl;
  }

  in_hashes_file.close();
  in_hashes_file.close();

  remove("tmp_hashes.bin");
  remove("tmp_counts.bin");

  for (uint i = 0; i < m.size(); i++) {
    REQUIRE(m[hashes[i]] == counts[i]);
  }

}

TEST_CASE("ParseKmerCountsAndCreateHDF5", "[HashTest]") {
  std::vector<uint> hashes;
  std::vector<uint> counts;
  btree::btree_map<uint, uint> hashed_counts;
  const string kmer1("AAAAA");
  const string kmer2("TTTTT");
  const uint kmer1_count = 84;
  const uint kmer2_count = 13;
  uint kmer1_pos = 0;
  uint kmer2_pos = 0;
  uint pos = 0, last;
  FQ::DataType type;
  string sample_hash_dataset_name("kmer_hash"), sample_counts_dataset_name("count");
  uint *hashes_arr, *counts_arr;

  std::vector<uint64_t> dims;
  param_struct params;

  params.k_val_start = 5;
  params.k_val_end = 5;
  params.hdf5_filename = "/home/darryl/Development/kmerge/tests/parse_example.h5";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/genome.test.fasta.gz";
  params.tmp_hashes_filename = "hashes.bin";
  params.tmp_counts_filename = "counts.bin";
  params.group_name = "/org1";
  params.hash_dataset_name =  "/org1/kmer_hash";
  params.counts_dataset_name = "/org1/count";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  params.seq_filename ="/home/darryl/Development/kmerge/tests/genome.test.fasta.gz";
  
  try {

    KMerge* kmerge = new KMerge(params.hdf5_filename, "lookup3", ".");
    params.kmerge = kmerge;

    bool success = kmerge->count_hashed_kmers_fasta(params, hashed_counts);
    REQUIRE(success == true);

    REQUIRE(hashed_counts.size() == 324);

    //make sure hashes are sorted and store in vectors
    for (btree::btree_map<uint, uint>::const_iterator m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
      if (m_iter != hashed_counts.begin()) {
        REQUIRE(m_iter->first > last);
      }
      hashes.push_back(m_iter->first);
      counts.push_back(m_iter->second);
      last = m_iter->first;
    }
    
    REQUIRE(hashes.size() == 324);
    REQUIRE(counts.size() == 324);

    for (std::vector<uint>::iterator v_iter = hashes.begin(); v_iter != hashes.end(); v_iter++) {
      REQUIRE(*v_iter >= 0);
      if (*v_iter == kmerge->hash_kmer(kmer1)) {
	kmer1_pos = pos;
      }
      if (*v_iter == kmerge->hash_kmer(kmer2)) {                         
        kmer2_pos = pos;            
      }  
      pos++;
    }
    REQUIRE(kmer1_pos == kmer2_pos);
    REQUIRE(hashes[kmer1_pos] == kmerge->hash_kmer(kmer2));
    REQUIRE(counts[kmer2_pos] == kmer1_count + kmer2_count);

    
    kmerge->add_dataset(params.hash_dataset_name, hashes.size(), &hashes[0], NULL);
    kmerge->add_dataset(params.counts_dataset_name, counts.size(), &counts[0], NULL);
    
    delete kmerge;

    FastQuery *sample = new FastQuery(params.hdf5_filename, FQ::FQ_HDF5);

    if (!(sample->getVariableInfo(sample_counts_dataset_name, sample_counts_dataset_name, dims, &type, params.group_name))) {
      cerr << "Cannot access sample variable information." << endl;
      exit(EXIT_FAILURE);
    }
    hashes_arr = new uint[dims[0]];
   
    if (!(sample->getData("kmer_hash", hashes_arr, &params.group_name[0]))) {
      cerr << "Cannot access sample hashes" << endl;
      exit(EXIT_FAILURE);
    }
    
    counts_arr = new uint[dims[0]];

    if(!(sample->getData("count", counts_arr, &params.group_name[0]))) {
      cerr << "Cannot access sample counts" << endl;
      exit(EXIT_FAILURE);
    }

    for (pos = 0; pos < dims[0]; pos++) {
      if (pos != 0) {
        REQUIRE(hashes_arr[pos] > last); //make sure sorting still intact
      }
      REQUIRE(hashes_arr[pos] == hashes[pos]);
      REQUIRE(counts_arr[pos] == counts[pos]);
      last = hashes_arr[pos];
    }
    
    
    delete hashes_arr;
    delete counts_arr;
    delete sample;

    /*if (remove(params.hdf5_filename.c_str()) != 0) {
      perror( "Error deleting file");
      }*/

  } catch (exception &e) {
    cout << e.what() << endl;
  } 
}


/*
TEST_CASE("ProblemGenomeTest", "[HashTest]") {
  param_struct params1, params2, params3;

  FQ::DataType type;
  std::vector<uint64_t> dims;
  std::string hash_dataset_name("kmer_hash"), counts_dataset_name("count");

  params1.k_val_start = 3;
  params1.k_val_end = 3;
  params1.hdf5_filename = "problem.h5";
  params1.seq_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/57623/57623.fasta.gz";
  params1.tmp_hashes_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/57623/hashes.bin";
  params1.tmp_counts_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/57623/counts.bin";
  params1.group_name = "/57623";
  params1.hash_dataset_name =  "/57623/kmer_hash";
  params1.counts_dataset_name = "/57623/count";
  params1.num_threads = (params1.k_val_end - params1.k_val_start) / 2 + 1;

  params2.k_val_start = 3;
  params2.k_val_end = 3;
  params2.hdf5_filename = "/home/darryl/Development/kmerge/tests/problem.h5";
  params2.seq_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/50385/50385.fasta.gz";
  params2.tmp_hashes_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/50385/hashes.bin";
  params2.tmp_counts_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/50385/counts.bin";
  params2.group_name = "/50385";
  params2.hash_dataset_name =  "/50385/kmer_hash";
  params2.counts_dataset_name = "/50385/count";
  params2.num_threads = (params2.k_val_end - params2.k_val_start) / 2 + 1;

  dlib::thread_pool tp(3);

  KMerge* kmerge = new KMerge(params1.hdf5_filename, "lookup3", ".");

  params1.kmerge = kmerge;
  params2.kmerge = kmerge;

  KMerge::BuilderTask t1(params1);
  KMerge::BuilderTask t2(params2);

  tp.add_task(t1, &KMerge::BuilderTask::execute);
  tp.add_task(t2, &KMerge::BuilderTask::execute);

  tp.wait_for_all_tasks();

  delete kmerge;

  HDF5 *sample = new HDF5(params1.hdf5_filename);

  if (!(sample->getVariableInfo(params1.counts_dataset_name, dims, &type))) {
    cerr << "Cannot access sample variable information."  << endl;
  }
  
  REQUIRE(dims[0] == MAX_UINT_VAL);
  dims.clear();

  if (!(sample->getVariableInfo(params2.counts_dataset_name, dims, &type))) {
    cerr << "Cannot access sample variable information."  << endl;
  }

  REQUIRE(dims[0] == MAX_UINT_VAL);


  delete sample;
   

  if (remove( params1.hdf5_filename.c_str() ) != 0) {
    perror( "Error deleting file");
  }

}
*/

TEST_CASE("ThreadedParseKmerCountsAndCreateHDF5FromFastq", "[HashTest]") {
  param_struct params;
  FQ::DataType type;
  std::vector<uint64_t> dims;

  params.k_val_start = 3;
  params.k_val_end = 7;
  params.hdf5_filename = "thread_fastq.h5";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz";
  params.tmp_hashes_filename = "";
  params.tmp_counts_filename = "";
  params.group_name = "/sample";
  params.hash_dataset_name =  "/sample/kmer_hash";
  params.counts_dataset_name = "/sample/count";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  KMerge *kmerge = new KMerge(params.hdf5_filename, "lookup3", ".");

  params.kmerge = kmerge;

  kmerge->build(params);

  delete kmerge;

  HDF5 *sample = new HDF5(params.hdf5_filename);

  if (!(sample->getVariableInfo(params.counts_dataset_name, dims, &type))) {
    cerr << "Cannot access sample variable information."  << endl;
  }
  
  REQUIRE(dims[0] == 8727);

  uint *hashes_arr = new uint[dims[0]];
   
  if (!(sample->getData(params.hash_dataset_name, hashes_arr))) {
    cerr << "Cannot access sample hashes" << endl;
    exit(EXIT_FAILURE);
  }
  
  uint *counts_arr = new uint[dims[0]];
  
  if(!(sample->getData(params.counts_dataset_name, counts_arr))) {
    cerr << "Cannot access sample counts" << endl;
    exit(EXIT_FAILURE);
  }

  dims.clear();

  delete [] counts_arr;
  delete [] hashes_arr;
  delete sample;
   

  if (remove( params.hdf5_filename.c_str() ) != 0) {
    perror( "Error deleting file");
  }
}

TEST_CASE("ThreadedParseKmerCountsAndCreateHDF5", "[HashTest]") {
  std::vector<uint> hashes, counts;
  uint *hashes_arr, *counts_arr, kmer1_count, kmer2_count, kmer1_pos, kmer2_pos, pos;
  param_struct params1, params2, params3;

  FQ::DataType type;
  std::vector<uint64_t> dims;
  std::string hash_dataset_name("kmer_hash"), counts_dataset_name("count"), line;
  ifstream in_file;

  const string kmer1("AAAAA");
  const string kmer2("GCGAT");

  params1.k_val_start = 5;
  params1.k_val_end = 5;
  params2.k_val_start = 5;
  params2.k_val_end = 5;
  params3.k_val_start = 5;
  params3.k_val_end = 5;


  params1.hdf5_filename = "thread_example.h5";
  params2.hdf5_filename = "thread_example.h5";
  params3.hdf5_filename = "thread_example.h5";

  params1.seq_filename = "/home/darryl/Development/kmerge/tests/208831/208831.fasta.gz";
  params1.tmp_hashes_filename = "/home/darryl/Development/kmerge/tests/208831/hashes.bin";
  params1.tmp_counts_filename = "/home/darryl/Development/kmerge/tests/208831/counts.bin";
  params1.group_name = "/208831";
  params1.hash_dataset_name =  "/208831/kmer_hash";
  params1.counts_dataset_name = "/208831/count";
  params1.num_threads = (params1.k_val_end - params1.k_val_start) / 2 + 1;


  params2.seq_filename = "/home/darryl/Development/kmerge/tests/209328/209328.fasta.gz";
  params2.tmp_hashes_filename = "/home/darryl/Development/kmerge/tests/209328/hashes.bin";
  params2.tmp_counts_filename = "/home/darryl/Development/kmerge/tests/209328/counts.bin";
  params2.group_name = "/209328";
  params2.hash_dataset_name = "/209328/kmer_hash";
  params2.counts_dataset_name = "/209328/count";
  params2.num_threads = (params2.k_val_end - params2.k_val_start) / 2 + 1;

  params3.seq_filename = "/home/darryl/Development/kmerge/tests/54095/54095.fasta.gz";
  params3.tmp_hashes_filename = "/home/darryl/Development/kmerge/tests/54095/hashes.bin";
  params3.tmp_counts_filename = "/home/darryl/Development/kmerge/tests/54095/counts.bin";
  params3.group_name = "/54095";
  params3.hash_dataset_name = "/54095/kmer_hash";
  params3.counts_dataset_name = "/54095/count";
  params3.num_threads = (params3.k_val_end - params3.k_val_start) / 2 + 1;

  uint thread_count = 3;

  dlib::thread_pool tp(thread_count);
  try {
    KMerge* kmerge = new KMerge(params1.hdf5_filename, "lookup3", ".");

    params1.kmerge = kmerge;
    params2.kmerge = kmerge;
    params3.kmerge = kmerge;

    KMerge::BuilderTask t1(params1);
    KMerge::BuilderTask t2(params2);
    KMerge::BuilderTask t3(params3);

    tp.add_task(t1, &KMerge::BuilderTask::execute);
    tp.add_task(t2, &KMerge::BuilderTask::execute);
    tp.add_task(t3, &KMerge::BuilderTask::execute);

    tp.wait_for_all_tasks();

    delete kmerge;

    FastQuery *sample = new FastQuery(params1.hdf5_filename, FQ::FQ_HDF5);
 
    if (!(sample->getVariableInfo(counts_dataset_name, counts_dataset_name, dims, &type, params1.group_name))) {
      cerr << "Cannot access sample variable information."  << endl;
    }

    //REQUIRE(dims[0] == 512);
    REQUIRE(dims[0] == MAX_UINT_VAL);

    /*    hashes_arr = new uint[dims[0]];

    if (!(sample->getData("kmer_hash", hashes_arr, &params1.group_name[0]))) {
      cerr << "Cannot access sample hashes" << endl;
    }
    */
    counts_arr = new uint[dims[0]];


    if(!(sample->getData("count", counts_arr, &params1.group_name[0]))) {
      cerr << "Cannot access sample counts" << endl;
      exit(EXIT_FAILURE);
    }

    kmer1_count = 6150 /*AAAAA*/ + 6021 /*TTTTT*/;
    kmer2_count = 10775 /*GCGAT*/ + 10855 /*ATCGC*/;
    /*
    for (pos = 0; pos < dims[0]; pos++) {
      if (hashes_arr[pos] == KMerge::hash_kmer(kmer1, LOOKUP3)) {
        kmer1_pos = pos;
      }
      if (hashes_arr[pos] == KMerge::hash_kmer(kmer2, LOOKUP3)) {
        kmer2_pos = pos;
      }
    }

    REQUIRE(hashes_arr[kmer1_pos] == KMerge::hash_kmer(kmer1, LOOKUP3));
    REQUIRE(hashes_arr[kmer2_pos] == KMerge::hash_kmer(kmer2, LOOKUP3));
    REQUIRE(counts_arr[kmer1_pos] == kmer1_count);
    REQUIRE(counts_arr[kmer2_pos] == kmer2_count);
    */
    //AAAAA/TTTTT : 1130773752
    REQUIRE(counts_arr[1130773752] == kmer1_count);
    //GCGAT/ATCGC : 1282497623
    REQUIRE(counts_arr[1282497623] == kmer2_count);

    delete [] counts_arr;
    //delete [] hashes_arr;

    in_file.open("/home/darryl/Development/kmerge/tests/208831/taxonomy.txt");

    while (std::getline(in_file, line)) {
      std::istringstream tokenizer(line);
      std::stringstream path;
      std::string classification;
      std::vector<std::string> vars;
      uint i = 0;
      while (!tokenizer.eof()) {
	std::string token;
	getline(tokenizer, token, '\t');
	if (i == 0) { // this is taxon                                                                                         
	  path << "/208831/taxonomy/" << token;
	} else { // this is classification                                                                                     
	  classification = token;
	}
	i++;
      }
      if(!sample->getAllVariables(vars, path.str())) {
	throw "Unable to get variables";
      }
      path << "/";
      vars[0].replace(0, path.str().length(), "");
      REQUIRE(classification == vars[0]);
    }

    in_file.close();
    /*
    hashes_arr = new uint[dims[0]];

    if (!(sample->getData("kmer_hash", hashes_arr, &params2.group_name[0]))) {
      cerr << "Cannot access sample hashes" << endl;
      exit(EXIT_FAILURE);
    }
    */
    counts_arr = new uint[dims[0]];


    if(!(sample->getData("count", counts_arr, &params2.group_name[0]))) {
      cerr << "Cannot access sample counts" << endl;
    }

    kmer1_count = 2147 /*AAAAA*/ + 1919 /*TTTTT*/;
    kmer2_count = 12082 /*GCGAT*/ + 12213 /*ATCGC*/;
    /*
    for (pos = 0; pos < dims[0]; pos++) {
      if (hashes_arr[pos] == KMerge::hash_kmer(kmer1, LOOKUP3)) {
        kmer1_pos = pos;
      }
      if (hashes_arr[pos] == KMerge::hash_kmer(kmer2, LOOKUP3)) {
        kmer2_pos = pos;
      }
    }
   
    REQUIRE(hashes_arr[kmer1_pos] == KMerge::hash_kmer(kmer1, LOOKUP3));
    REQUIRE(hashes_arr[kmer2_pos] == KMerge::hash_kmer(kmer2, LOOKUP3));
    REQUIRE(counts_arr[kmer1_pos] == kmer1_count);
    REQUIRE(counts_arr[kmer2_pos] == kmer2_count);
    */
    //AAAAA/TTTTT : 1130773752 
    REQUIRE(counts_arr[1130773752] == kmer1_count);
    //GCGAT/ATCGC : 1282497623
    REQUIRE(counts_arr[1282497623] == kmer2_count);
    
    delete [] counts_arr;
    //delete [] hashes_arr;

    in_file.open("/home/darryl/Development/kmerge/tests/209328/taxonomy.txt");

    while (std::getline(in_file, line)) {
      std::istringstream tokenizer(line);
      std::stringstream path;
      std::string classification;
      std::vector<std::string> vars;
      uint i = 0;
      while (!tokenizer.eof()) {
	std::string token;
        getline(tokenizer, token, '\t');
        if (i == 0) { // this is taxon                                                                                         
          path << "/209328/taxonomy/" << token;
        } else { // this is classification                                                                                     
          classification = token;
        }
        i++;
      }
      if(!sample->getAllVariables(vars, path.str())) {
        throw "Unable to get variables";
      }
      path << "/";
      vars[0].replace(0, path.str().length(), "");
      REQUIRE(classification == vars[0]);
    }
    
    in_file.close();
    /*
    hashes_arr = new uint[dims[0]];

    if (!(sample->getData("kmer_hash", hashes_arr, &params3.group_name[0]))) {
      cerr << "Cannot access sample hashes" << endl;
      exit(EXIT_FAILURE);
    }
    */
    counts_arr = new uint[dims[0]];


    if(!(sample->getData("count", counts_arr, &params3.group_name[0]))) {
      cerr << "Cannot access sample counts" << endl;
    }
    
    kmer1_count = 9896 /*AAAAA*/ + 9505 /*TTTTT*/;
    kmer2_count = 733 /*GCGAT*/ + 750 /*ATCGC*/;
    /*
    for (pos = 0; pos < dims[0]; pos++) {
      if (hashes_arr[pos] == KMerge::hash_kmer(kmer1, LOOKUP3)) {
	kmer1_pos = pos;
      }
      if (hashes_arr[pos] == KMerge::hash_kmer(kmer2, LOOKUP3)) {
        kmer2_pos = pos;
      }
    }

    REQUIRE(hashes_arr[kmer1_pos] == KMerge::hash_kmer(kmer1, LOOKUP3));
    REQUIRE(hashes_arr[kmer2_pos] == KMerge::hash_kmer(kmer2, LOOKUP3));
    REQUIRE(counts_arr[kmer1_pos] == kmer1_count);
    REQUIRE(counts_arr[kmer2_pos] == kmer2_count);
    */


    //AAAAA/TTTTT : 1130773752
    REQUIRE(counts_arr[1130773752] == kmer1_count);
    //GCGAT/ATCGC : 1282497623
    REQUIRE(counts_arr[1282497623] == kmer2_count);

    delete [] counts_arr;
    //delete [] hashes_arr;

    in_file.open("/home/darryl/Development/kmerge/tests/54095/taxonomy.txt");

    while (std::getline(in_file, line)) {
      std::istringstream tokenizer(line);
      std::stringstream path;
      std::string classification;
      std::vector<std::string> vars;
      uint i = 0;
      while (!tokenizer.eof()) {
	std::string token;
        getline(tokenizer, token, '\t');
        if (i == 0) { // this is taxon                                                                                         
          path << "/54095/taxonomy/" << token;
        } else { // this is classification                                                                                     
          classification = token;
	}
        i++;
      }
      if(!sample->getAllVariables(vars, path.str())) {
	throw "Unable to get variables";
      }
      path << "/";
      vars[0].replace(0, path.str().length(), "");
      REQUIRE(classification == vars[0]);
    }

    in_file.close();

    delete sample;

    
    if (remove( params1.hdf5_filename.c_str() ) != 0) {
      perror( "Error deleting file");
    }
  } catch (exception &e) {
    cout << e.what() << endl;
  }
    
}

TEST_CASE("TestHashingFunctions", "[HashTest]") {
  btree::btree_map<uint, uint> hashed_counts;
  btree::btree_map<uint, uint>::const_iterator m_iter;
  
  param_struct params;

  params.k_val_start = 5;
  params.k_val_end = 5;
  params.hdf5_filename = "hash.h5";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/genome.test.fasta.gz";
  params.tmp_hashes_filename = "hashes.bin";
  params.tmp_counts_filename = "counts.bin";
  params.group_name = "/org1";
  params.hash_dataset_name =  "/org1/kmer_hash";
  params.counts_dataset_name = "/org1/count";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  params.kmerge = new KMerge(params.hdf5_filename, "lookup3", ".");
  bool success = params.kmerge->count_hashed_kmers_fasta(params, hashed_counts);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->first <= (uint) MAX_UINT_VAL);
  }

  delete params.kmerge;

  params.kmerge = new KMerge(params.hdf5_filename, "spooky", ".");
  success = params.kmerge->count_hashed_kmers_fasta(params, hashed_counts);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->first <= (uint) MAX_UINT_VAL);
  }

  delete params.kmerge;

  params.kmerge = new KMerge(params.hdf5_filename, "city", ".");
  success = params.kmerge->count_hashed_kmers_fasta(params, hashed_counts);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->first <= (uint) MAX_UINT_VAL);
  }

  delete params.kmerge;

  params.kmerge = new KMerge(params.hdf5_filename, "murmur", ".");
  success = params.kmerge->count_hashed_kmers_fasta(params, hashed_counts);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->first <= (uint) MAX_UINT_VAL);
  }

  delete params.kmerge;


  if (remove("hash.h5" ) != 0) {
    perror( "Error deleting file");
  }
}

TEST_CASE("AddTaxonomyInfoToHDF5File", "[HDF5Test]") {
  std::string hdf5_filename("taxonomy.h5"), group("/54095"); 
  std::string path_root("/54095/taxonomy"), line;
  

  KMerge *kmerge = new KMerge(hdf5_filename, "lookup3", ".");

  kmerge->add_taxonomy(group);

  delete kmerge;

  HDF5 *test = new HDF5(hdf5_filename, false);

  ifstream in_file("/home/darryl/Development/kmerge/tests/54095/taxonomy.txt");

  while (std::getline(in_file, line)) {
    std::istringstream tokenizer(line);
    std::stringstream path;
    std::string classification;
    std::vector<std::string> vars;
    uint i = 0;
    while (!tokenizer.eof()) {
      std::string token;
      getline(tokenizer, token, '\t');
      if (i == 0) { // this is taxon                                                                                         
	path << path_root << "/" << token;
      } else { // this is classification                                                                                     
	classification = token;
      }
      i++;
    }
    if(!test->getAllVariables(path.str(), vars)) {
      throw "Unable to get variables";
    }
    path << "/";
    vars[0].replace(0, path.str().length(), "");
    REQUIRE(classification == vars[0]);
  }



  delete test;

  if (remove(&hdf5_filename[0]) != 0) {
    perror( "Error deleting file");
  }
}
