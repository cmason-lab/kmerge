#define CATCH_CONFIG_MAIN

#include <limits.h>
#include "H5Cpp.h"
#include "kmerge.h"
#include <math.h>
#include "queryProcessor.h"
#include "indexBuilder.h"
#include "fq.h"
#include "hdf5file.h"
#include "catch.hpp"
#include "armadillo"
#include <dlib/threads.h>
#include <seqan/sequence.h>
#include <seqan/alignment_free.h>
#include <libGkArrays/gkArrays.h>

using namespace seqan;
using namespace dlib;


TEST_CASE("CountKmersWithGKArraysTest", "[HashTest]") {
  std::map<std::string, uint> gka_counts;
  // Default k-mer length
  int k = 5;

  ifstream in_file ("genome.test.kmers.txt");
  std::string genome_file("genome.test.contig.fa.gz");

  // Building the index
  gkarrays::gkArrays *genome = new gkarrays::gkArrays(&genome_file[0], k, true, 0, true);
  
  // checking number of indexed sequences
  
  REQUIRE(genome->getNbTags() == 3);
  REQUIRE(genome->getGkCFALength() == 424);

  gkarrays::readIterator *read_iter = genome->getReads()->begin();
  uint i = 0, num_iterations = 0;
  bool all_kmers_found = false;
  while (!read_iter->isFinished() && !all_kmers_found) {
    uint *support = genome->getSupport(i);
    for(uint j = 0; j < genome->getSupportLength(i); j++) {
      num_iterations++;
      char *kmer = genome->getTagFactor(i, j, k);
      if (gka_counts.count(kmer) == 0) {
	gka_counts[kmer] = support[j];
      }
      if (genome->getGkCFALength() == gka_counts.size()) {
	cout << "Breaking early (iteration = " << num_iterations << ")" << endl;
	all_kmers_found = true;
	break;
      }
      delete [] kmer;
    }
    delete [] support;
    ++(*read_iter);
    i++;
  }

  // compare to jellyfish results
  // perl ../scripts/contiguous_fasta.pl <(zcat genome.test.fa.gz) | ~/bin/fastx_toolkit/bin/fasta_formatter -w 80 | gzip > genome.test.contig.fa.gz
  // ~/bin/jellyfish count -m 5 -o output -s 10000 -C <(zcat genome.test.contig.fa.gz)
  // ~/bin/jellyfish dump -c output > genome.test.kmers.txt

  std::string seq;
  uint count;
  if (in_file.is_open()) {
    while ( !in_file.eof() ) {
      in_file >> seq >> count;
      REQUIRE(gka_counts[seq] == count);
    }
    in_file.close();
  }

  /* Strategy:
   * Keep a map of the k-mers encountered from a genome, loop through each sequence in the genome
   * for every k-mer encountered check if we've seen it before, if we have, ignore it,
   * if we haven't, look at total times k-mer is seen (stored in "count" of "getTagsWithFactor")
   * record this number of times in k-mer counts for genome and put it in kmers_found map so that it is ignored next time
   * need to loop over all sequences in the file when doing this to make sure all k-mers are acconted for
   * STOP when the kmers_found map has the same number of elements as "getGkCFALength".
   */

  // Free memory.
  delete read_iter;
  delete genome;
}

TEST_CASE("ReverseComplementTest", "[HashTest]") {
  std::string test_seq("ACTAG");
  std::string rc_test_seq(test_seq.c_str());
  reverseComplement(rc_test_seq);
  REQUIRE(rc_test_seq == "CTAGT");
}

TEST_CASE("TestHashedKmersAndReverseComplementReturnSameHashVal", "[HashTest]") {
  std::string test_seq("ACTAG");
  std::string rc_test_seq(test_seq.c_str());
  reverseComplement(rc_test_seq);
  
  KMerge *kmerge = new KMerge("dummy.h5", "lookup3");
  
  uint hash1 = kmerge->hash_kmer(test_seq), hash2 = kmerge->hash_kmer(rc_test_seq);
  REQUIRE(hash1 == hash2);
  delete kmerge;
}


TEST_CASE("ParseKmerCountsAndCreateHDF5", "[HashTest]") {
  std::vector<uint> hashes;
  std::vector<uint> counts;
  std::map<uint, uint> hashed_counts;
  const string kmer1("AAAAA");
  const string kmer2("TTTTT");
  const uint kmer1_count = 84;
  const uint kmer2_count = 13;
  uint kmer1_pos = 0;
  uint kmer2_pos = 0;
  uint pos = 0, last;
  FQ::DataType type;
  const string sample_var_path("/org1");
  string sample_hash_dataset_name("kmer_hash"), sample_counts_dataset_name("count");
  uint *hashes_arr, *counts_arr;

  const std::string HDF5_FILE_NAME( "/home/darryl/Development/kmerge/tests/parse_example.h5" );
  const std::string HASH_DATASET_NAME( "/org1/kmer_hash" );
  const std::string COUNT_DATASET_NAME( "/org1/count" );
  std::vector<uint64_t> dims;
  uint k = 5;
  std::string filename("genome.test.contig.fa.gz");

  try {

    KMerge* kmerge = new KMerge(HDF5_FILE_NAME, "lookup3");

    bool success = kmerge->count_hashed_kmers(filename, k, hashed_counts);
    REQUIRE(success == true);

    REQUIRE(hashed_counts.size() == 324);

    //make sure hashes are sorted and store in vectors
    for (std::map<uint, uint>::const_iterator m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
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

    
    kmerge->add_dataset(HASH_DATASET_NAME, hashes.size(), &hashes[0]);
    kmerge->add_dataset(COUNT_DATASET_NAME, counts.size(), &counts[0]);
    
    delete kmerge;

    FastQuery *sample = new FastQuery(HDF5_FILE_NAME, FQ::FQ_HDF5);

    if (!(sample->getVariableInfo(sample_hash_dataset_name, sample_hash_dataset_name, dims, &type, sample_var_path))) {
      cerr << "Cannot access sample variable information." << endl;
      exit(EXIT_FAILURE);
    }
    hashes_arr = new uint[dims[0]];
   
    if (!(sample->getData("kmer_hash", hashes_arr, &sample_var_path[0]))) {
      cerr << "Cannot access sample hashes" << endl;
      exit(EXIT_FAILURE);
    }
    
    counts_arr = new uint[dims[0]];

    if(!(sample->getData("count", counts_arr, &sample_var_path[0]))) {
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

    if (remove("/home/darryl/Development/kmerge/tests/parse_example.h5" ) != 0) {
      perror( "Error deleting file");
    }

  } catch (exception &e) {
    cout << e.what() << endl;
  } 
}


TEST_CASE("ThreadedParseKmerCountsAndCreateHDF5", "[HashTest]") {
  std::vector<uint> hashes, counts;
  uint *hashes_arr, *counts_arr, kmer1_count, kmer2_count, kmer1_pos, kmer2_pos, pos;
  param_struct params1, params2, params3;

  FQ::DataType type;
  std::vector<uint64_t> dims;
  std::string hash_dataset_name("kmer_hash"), counts_dataset_name("count");

  const string kmer1("AAAAA");
  const string kmer2("GCGAT");

  params1.k_val_start = 5;
  params1.k_val_end = 5;
  params2.k_val_start = 5;
  params2.k_val_end = 5;
  params3.k_val_start = 5;
  params3.k_val_end = 5;


  params1.hdf5_filename = "/home/darryl/Development/kmerge/tests/thread_example.h5";
  params2.hdf5_filename = "/home/darryl/Development/kmerge/tests/thread_example.h5";
  params3.hdf5_filename = "/home/darryl/Development/kmerge/tests/thread_example.h5";

  params1.seq_filename = "/home/darryl/Development/kmerge/tests/sandbox/15667/15667.fasta.gz";
  params1.group_name = "/15667";
  params1.hash_dataset_name =  "/15667/kmer_hash";
  params1.counts_dataset_name = "/15667/count";

  params2.seq_filename = "/home/darryl/Development/kmerge/tests/sandbox/165193/165193.fasta.gz";
  params2.group_name = "/165193";
  params2.hash_dataset_name = "/165193/kmer_hash";
  params2.counts_dataset_name = "/165193/count";

  params3.seq_filename = "/home/darryl/Development/kmerge/tests/sandbox/29363/29363.fasta.gz";
  params3.group_name = "/29363";
  params3.hash_dataset_name = "/29363/kmer_hash";
  params3.counts_dataset_name = "/29363/count";

  uint thread_count = 3;

  thread_pool tp(thread_count);
  try {
    KMerge* kmerge = new KMerge(params1.hdf5_filename, "lookup3");

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


    FastQuery *sample = new FastQuery(params1.hdf5_filename, FQ::FQ_HDF5);

    if (!(sample->getVariableInfo(hash_dataset_name, hash_dataset_name, dims, &type, params1.group_name))) {
      cerr << "Cannot access sample variable information." << endl;
      exit(EXIT_FAILURE);
    }

    REQUIRE(dims[0] == 512);

    hashes_arr = new uint[dims[0]];

    if (!(sample->getData("kmer_hash", hashes_arr, &params1.group_name[0]))) {
      cerr << "Cannot access sample hashes" << endl;
      exit(EXIT_FAILURE);
    }

    counts_arr = new uint[dims[0]];


    if(!(sample->getData("count", counts_arr, &params1.group_name[0]))) {
      cerr << "Cannot access sample counts" << endl;
      exit(EXIT_FAILURE);
    }

    kmer1_count = 10816 /*AAAAA*/ + 11439 /*TTTTT*/;
    kmer2_count = 7617 /*GCGAT*/ + 7672 /*ATCGC*/;

    for (pos = 0; pos < dims[0]; pos++) {
      if (hashes_arr[pos] == kmerge->hash_kmer(kmer1)) {
        kmer1_pos = pos;
      }
      if (hashes_arr[pos] == kmerge->hash_kmer(kmer2)) {
        kmer2_pos = pos;
      }
    }

    REQUIRE(hashes_arr[kmer1_pos] == kmerge->hash_kmer(kmer1));
    REQUIRE(hashes_arr[kmer2_pos] == kmerge->hash_kmer(kmer2));
    REQUIRE(counts_arr[kmer1_pos] == kmer1_count);
    REQUIRE(counts_arr[kmer2_pos] == kmer2_count);


    delete [] counts_arr;
    delete [] hashes_arr;

    hashes_arr = new uint[dims[0]];

    if (!(sample->getData("kmer_hash", hashes_arr, &params2.group_name[0]))) {
      cerr << "Cannot access sample hashes" << endl;
      exit(EXIT_FAILURE);
    }

    counts_arr = new uint[dims[0]];


    if(!(sample->getData("count", counts_arr, &params2.group_name[0]))) {
      cerr << "Cannot access sample counts" << endl;
    }

    kmer1_count = 6910 /*AAAAA*/ + 6830 /*TTTTT*/;
    kmer2_count = 1881 /*GCGAT*/ + 1926 /*ATCGC*/;

    for (pos = 0; pos < dims[0]; pos++) {
      if (hashes_arr[pos] == kmerge->hash_kmer(kmer1)) {
        kmer1_pos = pos;
      }
      if (hashes_arr[pos] == kmerge->hash_kmer(kmer2)) {
        kmer2_pos = pos;
      }
    }

    REQUIRE(hashes_arr[kmer1_pos] == kmerge->hash_kmer(kmer1));
    REQUIRE(hashes_arr[kmer2_pos] == kmerge->hash_kmer(kmer2));
    REQUIRE(counts_arr[kmer1_pos] == kmer1_count);
    REQUIRE(counts_arr[kmer2_pos] == kmer2_count);

    delete [] counts_arr;
    delete [] hashes_arr;

    hashes_arr = new uint[dims[0]];

    if (!(sample->getData("kmer_hash", hashes_arr, &params3.group_name[0]))) {
      cerr << "Cannot access sample hashes" << endl;
      exit(EXIT_FAILURE);
    }

    counts_arr = new uint[dims[0]];


    if(!(sample->getData("count", counts_arr, &params3.group_name[0]))) {
      cerr << "Cannot access sample counts" << endl;
    }

    kmer1_count = 11579 /*AAAAA*/ + 12074 /*TTTTT*/;
    kmer2_count = 234 /*GCGAT*/ + 234 /*ATCGC*/;

    for (pos = 0; pos < dims[0]; pos++) {
      if (hashes_arr[pos] == kmerge->hash_kmer(kmer1)) {
	kmer1_pos = pos;
      }
      if (hashes_arr[pos] == kmerge->hash_kmer(kmer2)) {
        kmer2_pos = pos;
      }
    }

    REQUIRE(hashes_arr[kmer1_pos] == kmerge->hash_kmer(kmer1));
    REQUIRE(hashes_arr[kmer2_pos] == kmerge->hash_kmer(kmer2));
    REQUIRE(counts_arr[kmer1_pos] == kmer1_count);
    REQUIRE(counts_arr[kmer2_pos] == kmer2_count);

    delete [] counts_arr;
    delete [] hashes_arr;

    delete sample;

    delete kmerge;
    
    if (remove("/home/darryl/Development/kmerge/tests/thread_example.h5" ) != 0) {
      perror( "Error deleting file");
    }
  } catch (exception &e) {
    cout << e.what() << endl;
  }
    
}

TEST_CASE("TestHashingFunctions", "[HashTest]") {
  std::string filename("genome.test.contig.fa.gz"), hdf5_filename("/home/darryl/Development/kmerge/tests/hash.h5");
  KMerge* kmerge;
  uint k=5;
  std::map<uint, uint> hashed_counts;
  std::map<uint, uint>::const_iterator m_iter;

  kmerge = new KMerge(hdf5_filename, "lookup3");
  bool success = kmerge->count_hashed_kmers(filename, k, hashed_counts);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->first <= (uint) MAX_UINT_VAL);
  }

  delete kmerge;

  kmerge = new KMerge(hdf5_filename, "spooky");
  success = kmerge->count_hashed_kmers(filename, k, hashed_counts);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->first <= (uint) MAX_UINT_VAL);
  }

  delete kmerge;

  kmerge = new KMerge(hdf5_filename, "city");
  success = kmerge->count_hashed_kmers(filename, k, hashed_counts);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->first <= (uint) MAX_UINT_VAL);
  }

  delete kmerge;

  kmerge = new KMerge(hdf5_filename, "murmur");
  success = kmerge->count_hashed_kmers(filename, k, hashed_counts);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->first <= (uint) MAX_UINT_VAL);
  }


  if (remove("/home/darryl/Development/kmerge/tests/hash.h5" ) != 0) {
    perror( "Error deleting file");
  }
}
