#define CATCH_CONFIG_MAIN

#include <limits.h>
#include "kmerge.h"
#include <math.h>
#include "fq.h"
#include "hdf5file.h"
#include "catch.hpp"
#include <sstream>
#include <dlib/threads.h>
#include <dlib/serialize.h>
#include <dlib/logger.h>
#include "cpp-btree/btree_map.h"
#include <iomanip>


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

TEST_CASE("CountdKmersTest", "[HashTest]") {

  int k = 5, l;
  std::map<std::string, uint> kmer_counts, jf_counts;
  kseq_t *seq;
  gzFile fp;

  fp = gzopen("genome.test.fasta.gz", "r");
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

TEST_CASE("CountKmersInFastqFile", "[HashTest]") {

  int k = 7, l;
  std::map<std::string, uint> kmer_counts, jf_counts;
  kseq_t *seq;
  gzFile fp;

  fp = gzopen("sample/sample.fastq.gz", "r");
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

  ifstream in_file ("sample/sample.k7.txt");

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
  params.hdf5_filename = "/home/darryl/Development/kmerge/tests/fastq.h5";
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
  params.hdf5_filename = "/home/darryl/Development/kmerge/tests/problem.h5";
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
  params.seq_filename = "genome.test.fasta.gz";
  params.tmp_hashes_filename = "hashes.bin";
  params.tmp_counts_filename = "counts.bin";
  params.group_name = "/org1";
  params.hash_dataset_name =  "/org1/kmer_hash";
  params.counts_dataset_name = "/org1/count";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  params.seq_filename ="genome.test.fasta.gz";
  
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

    if (!(sample->getVariableInfo(sample_hash_dataset_name, sample_hash_dataset_name, dims, &type, params.group_name))) {
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

    if (remove(params.hdf5_filename.c_str()) != 0) {
      perror( "Error deleting file");
    }

  } catch (exception &e) {
    cout << e.what() << endl;
  } 
}

TEST_CASE("ProblemGenomeTest", "[HashTest]") {
  param_struct params1, params2, params3;

  FQ::DataType type;
  std::vector<uint64_t> dims;
  std::string hash_dataset_name("kmer_hash"), counts_dataset_name("count");

  params1.k_val_start = 3;
  params1.k_val_end = 3;
  params1.hdf5_filename = "/home/darryl/Development/kmerge/tests/problem.h5";
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
  /*
  params3.k_val_start = 3;
  params3.k_val_end = 3;
  params3.hdf5_filename = "/home/darryl/Development/kmerge/tests/problem.h5";
  params3.seq_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/202751/202751.fasta.gz";
  params3.tmp_hashes_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/202751/hashes.bin";
  params3.tmp_counts_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/202751/counts.bin";
  params3.group_name = "/202751";
  params3.hash_dataset_name =  "/202751/kmer_hash";
  params3.counts_dataset_name = "/202751/count";
  params3.num_threads = (params3.k_val_end - params3.k_val_start) / 2 + 1;
  */
  dlib::thread_pool tp(3);

  KMerge* kmerge = new KMerge(params1.hdf5_filename, "lookup3", ".");

  params1.kmerge = kmerge;
  params2.kmerge = kmerge;
  //params3.kmerge = kmerge;

  KMerge::BuilderTask t1(params1);
  KMerge::BuilderTask t2(params2);
  //KMerge::BuilderTask t3(params3);

  tp.add_task(t1, &KMerge::BuilderTask::execute);
  tp.add_task(t2, &KMerge::BuilderTask::execute);
  //tp.add_task(t3, &KMerge::BuilderTask::execute);

  tp.wait_for_all_tasks();

  delete kmerge;

  HDF5 *sample = new HDF5(params1.hdf5_filename);

  if (!(sample->getVariableInfo(params1.hash_dataset_name, dims, &type))) {
    cerr << "Cannot access sample variable information."  << endl;
  }
  
  REQUIRE(dims[0] == 32);
  dims.clear();

  if (!(sample->getVariableInfo(params2.hash_dataset_name, dims, &type))) {
    cerr << "Cannot access sample variable information."  << endl;
  }

  REQUIRE(dims[0] == 32);


  delete sample;
   

  if (remove( params1.hdf5_filename.c_str() ) != 0) {
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


  params1.hdf5_filename = "/home/darryl/Development/kmerge/tests/thread_example.h5";
  params2.hdf5_filename = "/home/darryl/Development/kmerge/tests/thread_example.h5";
  params3.hdf5_filename = "/home/darryl/Development/kmerge/tests/thread_example.h5";

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
  params.hdf5_filename = "/home/darryl/Development/kmerge/tests/hash.h5";
  params.seq_filename = "genome.test.fasta.gz";
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
