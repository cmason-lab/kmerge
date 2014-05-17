#define CATCH_CONFIG_MAIN

#include <limits.h>
#include "kmerge.h"
#include <math.h>
#include "fq.h"
#include "hdf5file.h"
#include "catch.hpp"
#include <sstream>
#include <dlib/threads.h>
#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;
using namespace dlib;


TEST_CASE("CountHashedKmers", "[HashTest") {
  int k = 5;
  std::map<std::string, uint> kmer_counts, jf_counts;
  seqan::SequenceStream seq_io("genome.test.fasta.gz");
  REQUIRE(isGood(seq_io) == true);
  // Read file record-wise.
  seqan::CharString id;
  seqan::Dna5String seq;
  std::stringstream kmer;
  while (!atEnd(seq_io)) {
    REQUIRE(readRecord(id, seq, seq_io) == 0);
    for (int i = 0; i < length(seq) - k + 1; i++) {
      seqan::Infix<seqan::Dna5String>::Type sub_seq(seq, i, i+k);
      
      kmer.str("");
      kmer << sub_seq;

      if(kmer.str().find("N") != std::string::npos) { // skip kmers containing Ns
	continue;
      }
      if (kmer_counts.find(kmer.str()) != kmer_counts.end()) {
	kmer_counts[kmer.str()]++;
      } else {
	kmer_counts[kmer.str()] = 1;
      }
    }
  }


  // compare to jellyfish results
  // perl ../scripts/contiguous_fasta.pl <(zcat genome.test.fa.gz) | ~/bin/fastx_toolkit/bin/fasta_formatter -w 80 | gzip > genome.test.contig.fa.gz                                            // ~/bin/jellyfish count -m 5 -o output -s 10000 <(zcat genome.test.contig.fa.gz)                                                  
  // ~/bin/jellyfish dump -c output > genome.test.kmers.txt                                                                          

  ifstream in_file ("genome.test.kmers.txt");

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
  
  KMerge *kmerge = new KMerge("dummy.h5", "lookup3", ".");
  
  uint hash1 = kmerge->hash_kmer(test_seq), hash2 = kmerge->hash_kmer(rc_test_seq);
  REQUIRE(hash1 == hash2);
  delete kmerge;
}

TEST_CASE("SerializeAndDeserializeMap", "[SerializeTest]") {
  std::map<uint, uint> m;
  ofstream out_file("test_serialize.txt");
  uint hash, count;
  std::vector<uint> hashes, counts;

  m[0] = 4;
  m[1] = 1;
  m[2] = 0;
  m[3] = 2;

  for (std::map<uint, uint>::const_iterator m_iter = m.begin(); m_iter != m.end(); m_iter++) {
    out_file << m_iter->first << "\t" << m_iter->second << endl;
  }
  out_file.close();

  ifstream in_file("test_serialize.txt");

  while (in_file >> hash >> count) {
    hashes.push_back(hash);
    counts.push_back(count);
  }
  in_file.close();
  remove("test_serialize.txt");

  for (uint i = 0; i < m.size(); i++) {
    REQUIRE(m[hashes[i]] == counts[i]);
  }
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
  std::string filename("genome.test.fasta.gz");

  try {

    KMerge* kmerge = new KMerge(HDF5_FILE_NAME, "lookup3", ".");

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


  params1.hdf5_filename = "/home/darryl/Development/kmerge/tests/thread_example.h5";
  params2.hdf5_filename = "/home/darryl/Development/kmerge/tests/thread_example.h5";
  params3.hdf5_filename = "/home/darryl/Development/kmerge/tests/thread_example.h5";

  params1.seq_filename = "/home/darryl/Development/kmerge/tests/208831/208831.fasta.gz";
  params1.group_name = "/208831";
  params1.tmp_filename = "/home/darryl/Development/kmerge/tests/208831/tmp.txt";
  params1.hash_dataset_name =  "/208831/kmer_hash";
  params1.counts_dataset_name = "/208831/count";

  params2.seq_filename = "/home/darryl/Development/kmerge/tests/209328/209328.fasta.gz";
  params2.group_name = "/209328";
  params2.tmp_filename = "/home/darryl/Development/kmerge/tests/209328/tmp.txt";
  params2.hash_dataset_name = "/209328/kmer_hash";
  params2.counts_dataset_name = "/209328/count";

  params3.seq_filename = "/home/darryl/Development/kmerge/tests/54095/54095.fasta.gz";
  params3.group_name = "/54095";
  params3.tmp_filename = "/home/darryl/Development/kmerge/tests/54095/tmp.txt";
  params3.hash_dataset_name = "/54095/kmer_hash";
  params3.counts_dataset_name = "/54095/count";

  uint thread_count = 3;

  thread_pool tp(thread_count);
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

    if (!(sample->getVariableInfo(hash_dataset_name, hash_dataset_name, dims, &type, params1.group_name))) {
      cerr << "Cannot access sample variable information."  << endl;
    }

    REQUIRE(dims[0] == 512);

    hashes_arr = new uint[dims[0]];

    if (!(sample->getData("kmer_hash", hashes_arr, &params1.group_name[0]))) {
      cerr << "Cannot access sample hashes" << endl;
    }

    counts_arr = new uint[dims[0]];


    if(!(sample->getData("count", counts_arr, &params1.group_name[0]))) {
      cerr << "Cannot access sample counts" << endl;
      exit(EXIT_FAILURE);
    }

    kmer1_count = 6150 /*AAAAA*/ + 6021 /*TTTTT*/;
    kmer2_count = 10775 /*GCGAT*/ + 10855 /*ATCGC*/;

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


    delete [] counts_arr;
    delete [] hashes_arr;

    in_file.open("./208831/taxonomy.txt");

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

    hashes_arr = new uint[dims[0]];

    if (!(sample->getData("kmer_hash", hashes_arr, &params2.group_name[0]))) {
      cerr << "Cannot access sample hashes" << endl;
      exit(EXIT_FAILURE);
    }

    counts_arr = new uint[dims[0]];


    if(!(sample->getData("count", counts_arr, &params2.group_name[0]))) {
      cerr << "Cannot access sample counts" << endl;
    }

    kmer1_count = 2147 /*AAAAA*/ + 1919 /*TTTTT*/;
    kmer2_count = 12082 /*GCGAT*/ + 12213 /*ATCGC*/;

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

    delete [] counts_arr;
    delete [] hashes_arr;

    in_file.open("./209328/taxonomy.txt");

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

    hashes_arr = new uint[dims[0]];

    if (!(sample->getData("kmer_hash", hashes_arr, &params3.group_name[0]))) {
      cerr << "Cannot access sample hashes" << endl;
      exit(EXIT_FAILURE);
    }

    counts_arr = new uint[dims[0]];


    if(!(sample->getData("count", counts_arr, &params3.group_name[0]))) {
      cerr << "Cannot access sample counts" << endl;
    }
    
    kmer1_count = 9896 /*AAAAA*/ + 9505 /*TTTTT*/;
    kmer2_count = 733 /*GCGAT*/ + 750 /*ATCGC*/;
  
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

    delete [] counts_arr;
    delete [] hashes_arr;

    in_file.open("./54095/taxonomy.txt");

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

  kmerge = new KMerge(hdf5_filename, "lookup3", ".");
  bool success = kmerge->count_hashed_kmers(filename, k, hashed_counts);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->first <= (uint) MAX_UINT_VAL);
  }

  delete kmerge;

  kmerge = new KMerge(hdf5_filename, "spooky", ".");
  success = kmerge->count_hashed_kmers(filename, k, hashed_counts);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->first <= (uint) MAX_UINT_VAL);
  }

  delete kmerge;

  kmerge = new KMerge(hdf5_filename, "city", ".");
  success = kmerge->count_hashed_kmers(filename, k, hashed_counts);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->first <= (uint) MAX_UINT_VAL);
  }

  delete kmerge;

  kmerge = new KMerge(hdf5_filename, "murmur", ".");
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

TEST_CASE("AddTaxonomyInfoToHDF5File", "[HDF5Test]") {
  std::string hdf5_filename("/home/darryl/Development/kmerge/tests/taxonomy.h5"), group("/54095"); 
  std::string path_root("/54095/taxonomy"), line;
  

  KMerge *kmerge = new KMerge(hdf5_filename, "lookup3", ".");

  kmerge->add_taxonomy(group);

  delete kmerge;

  HDF5 *test = new HDF5(hdf5_filename, false);

  ifstream in_file("./54095/taxonomy.txt");

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
