#define CATCH_CONFIG_MAIN

#include <limits.h>
#include "kmerge.h"
#include <math.h>
#include "catch.hpp"
#include <sstream>
#include <dlib/threads.h>
#include <dlib/serialize.h>
#include <dlib/logger.h>
#include <iomanip>

TEST_CASE("CompressHashesTest", "CompressionTest") {
  size_t N = 10 * 1000;
  std::vector<uint32_t> mydata(N);
  for(uint32_t i = 0; i < N;i += 150) mydata[i] = i;

  std::vector<uint32_t> compressed_output = KMerge::compress(mydata);

  std::cout<<std::setprecision(3);
  std::cout<<"You are using " << 32.0 * static_cast<double>(compressed_output.size()) /
    static_cast<double>(mydata.size()) <<" bits per integer. "<<std::endl;
  

  std::vector<uint32_t> mydataback = KMerge::uncompress(compressed_output, N);

  for (uint i = 0; i <= mydata.size(); i++) {
    REQUIRE(mydata[i] == mydataback[i]);
  }

}


TEST_CASE("LevelDBBasicsTest", "HashStorageTest") {
  std::string key1("test"), value;
  uint value_out = 1, value_in;
  std::stringstream ss_in, ss_out;
  dlib::serialize(value_out, ss_out);
  leveldb::DB* db;
  leveldb::Options options;
  options.create_if_missing = true;
  leveldb::Status status = leveldb::DB::Open(options, "./testdb", &db);
  REQUIRE(status.ok() == true);

  leveldb::Status s = db->Put(leveldb::WriteOptions(), key1, ss_out.str());
  REQUIRE(s.ok() == true);
  s = db->Get(leveldb::ReadOptions(), key1, &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(value_in, ss_in);
  REQUIRE(value_in == 1);
  delete db;

  system("rm -r ./testdb");
}

TEST_CASE("LevelDBAndFastPForTest", "HashStorageTest") {
  std::string key1("test"), value;
  std::stringstream ss_in, ss_out;
  leveldb::DB* db;
  leveldb::Options options;
  options.create_if_missing = true;

  size_t N = 10 * 1000;
  std::vector<uint32_t> mydata(N);
  for(uint32_t i = 0; i < N;i += 150) mydata[i] = i;

  std::vector<uint32_t> compressed_output = KMerge::compress(mydata), data_in;

  dlib::serialize(compressed_output, ss_out);
  
  leveldb::Status status = leveldb::DB::Open(options, "./testdb", &db);
  REQUIRE(status.ok() == true);

  leveldb::Status s = db->Put(leveldb::WriteOptions(), key1, ss_out.str());
  REQUIRE(s.ok() == true);


  s = db->Get(leveldb::ReadOptions(), key1, &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(data_in, ss_in);
  

  std::vector<uint32_t> mydataback = KMerge::uncompress(compressed_output, N);

  for (uint i = 0; i <= mydata.size(); i++) {
    REQUIRE(mydata[i] == mydataback[i]);
  }

  delete db;

  system("rm -r ./testdb");

}

TEST_CASE("ChainedHashMap", "ConcurrentTest") {
  ulib::chain_hash_map<uint, uint> c_map(100000000);

  c_map[0] = 0;
  c_map[1] = 1;
  c_map[2] = 2;
  c_map[3] = 3;

  uint i = 0;
  for (ulib::chain_hash_map<uint,uint>::const_iterator iter = c_map.begin(); iter != c_map.end(); iter++) {
    REQUIRE(c_map[i] == i);
    i++;
  }

  c_map[0] += c_map[3];

  REQUIRE(c_map[0] == c_map[3]);

  REQUIRE(c_map[10] == 0);

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

TEST_CASE("CountKmersInFastqFileBasic", "[HashTest]") {

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
  std::map<uint, uint> jf_hashed_counts;
  ulib::chain_hash_map<uint,uint> hashed_counts(100000000);
  param_struct params;
  std::vector<std::string> files;

  params.k_val_start = 3;
  params.k_val_end = 7;
  params.db_filename = "/home/darryl/Development/kmerge/tests/fasta.db";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/208831/208831.fasta.gz";
  params.group_name = "208831";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;  


  //~/bin/kanalyze-0.9.5/count -d 20 -f fastagz -k 3 -o 208831/sample.k3.txt 208831/208831.fasta.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastagz -k 5 -o 208831/sample.k5.txt 208831/208831.fasta.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastagz -k 7 -o 208831/sample.k7.txt 208831/208831.fasta.gz

  files.push_back("/home/darryl/Development/kmerge/tests/208831/sample.k3.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/208831/sample.k5.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/208831/sample.k7.txt");

  KMerge *kmerge = new KMerge(params.db_filename.c_str(), "lookup3", ".");

  params.kmerge = kmerge;


  bool success = kmerge->count_hashed_kmers(params, hashed_counts, true);
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
  
  for (std::map<uint, uint>::const_iterator b_it = jf_hashed_counts.begin(); b_it != jf_hashed_counts.end(); b_it++) {
    REQUIRE(hashed_counts[b_it->first] == b_it->second);
  }

}

TEST_CASE("CountHashedKmersInFastqFile", "[HashTest]") {
  ulib::chain_hash_map<uint,uint> hashed_counts(100000000);
  btree::btree_map<uint, uint> jf_hashed_counts;
  param_struct params;
  std::vector<std::string> files;

  params.k_val_start = 3;
  params.k_val_end = 7;
  params.db_filename = "fastq.db";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz";
  params.group_name = "sample";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;  

  //~/bin/kanalyze-0.9.5/count -d 20 -f fastqgz -k 3 -o sample/sample.k3.txt sample/sample.fastq.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastqgz -k 5 -o sample/sample.k5.txt sample/sample.fastq.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastqgz -k 7 -o sample/sample.k7.txt sample/sample.fastq.gz

  files.push_back("/home/darryl/Development/kmerge/tests/sample/sample.k3.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/sample/sample.k5.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/sample/sample.k7.txt");

  KMerge *kmerge = new KMerge(params.db_filename.c_str(), "lookup3", ".");

  params.kmerge = kmerge;


  bool success = kmerge->count_hashed_kmers(params, hashed_counts, false);
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

}


TEST_CASE("CountHashedKmersInParallelFasta", "[HashTest]") {
  ulib::chain_hash_map<uint, uint> hashed_counts(100000000);
  int l;
  param_struct params;
  
  params.k_val_start = 3;
  params.k_val_end = 11;
  params.db_filename = "dummy.db";
  params.seq_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/208831/208831.fasta.gz";
  params.group_name = "208831";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
  
  KMerge *kmerge = new KMerge(params.db_filename, "lookup3", ".");
  params.kmerge = kmerge;

  KMerge::CountAndHashSeq func(params, hashed_counts, true);
  dlib::parallel_for(params.num_threads, (params.k_val_start + 1) / 2, (params.k_val_end + 1) / 2 + 1, func);
    
  delete kmerge;
  
  REQUIRE(hashed_counts.size() == 1591922);
  

}

TEST_CASE("CountHashedKmersInParallelFastq", "[HashTest]") {
  ulib::chain_hash_map<uint, uint> hashed_counts(100000000);
  KMerge *kmerge = new KMerge("dummy.db", "lookup3", ".");
  param_struct params;

  params.k_val_start = 3;
  params.k_val_end = 7;
  params.db_filename = "dummy.db";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz";
  params.group_name = "sample";
  params.kmerge = kmerge;
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;


  KMerge::CountAndHashSeq func(params, hashed_counts, false);
  dlib::parallel_for(params.num_threads, (params.k_val_start + 1) / 2, (params.k_val_end + 1) / 2 + 1, func);

  delete kmerge;

  REQUIRE(hashed_counts.size() == 8727);
}

TEST_CASE("TestHashedKmersAndReverseComplementReturnSameHashVal", "[HashTest]") {
  std::string test_seq("ACTAG");
  std::string rc_test_seq = KMerge::rev_comp(test_seq);
  
  REQUIRE(rc_test_seq == "CTAGT");
  KMerge *kmerge = new KMerge("dummy.db", "lookup3", ".");
  
  uint hash1 = kmerge->hash_kmer(test_seq), hash2 = kmerge->hash_kmer(rc_test_seq);
  REQUIRE(hash1 == hash2);
  delete kmerge;
}

TEST_CASE("ParseKmerCountsAndCreateDB", "[HashTest]") {
  std::vector<uint> hashes_out, comp_hashes;
  std::vector<uint> counts_out, comp_counts;
  ulib::chain_hash_map<uint, uint> hashed_counts(1000000);
  const string kmer1("AAAAA");
  const string kmer2("TTTTT");
  const uint kmer1_count = 84;
  const uint kmer2_count = 13;
  uint kmer1_pos = 0;
  uint kmer2_pos = 0;
  uint pos = 0;
  std::string value;
  param_struct params;
  std::stringstream ss_in, ss_delete;


  params.k_val_start = 5;
  params.k_val_end = 5;
  params.db_filename = "/home/darryl/Development/kmerge/tests/parse_example.db";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/genome.test.fasta.gz";
  params.group_name = "org1";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
  params.priority = 1;
  params.lock_filename = params.db_filename + std::string(".lck");

  params.seq_filename ="/home/darryl/Development/kmerge/tests/genome.test.fasta.gz";
  

  KMerge* kmerge = new KMerge(params.db_filename, "lookup3", ".");
  params.kmerge = kmerge;
  
  bool success = kmerge->count_hashed_kmers(params, hashed_counts, true);
  REQUIRE(success == true);
  REQUIRE(hashed_counts.size() == 324);

  //store hashes in vectors
  for (ulib::chain_hash_map<uint, uint>::const_iterator m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    hashes_out.push_back(m_iter.key());
    counts_out.push_back(m_iter.value());
  }
    
  REQUIRE(hashes_out.size() == 324);
  REQUIRE(counts_out.size() == 324);

  for (std::vector<uint>::iterator v_iter = hashes_out.begin(); v_iter != hashes_out.end(); v_iter++) {
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
  REQUIRE(hashes_out[kmer1_pos] == kmerge->hash_kmer(kmer2));
  REQUIRE(counts_out[kmer2_pos] == kmer1_count + kmer2_count);

  kmerge->add_dataset_size(hashes_out.size(), params.group_name + std::string("|size"));  
  kmerge->add_dataset(hashes_out, params.group_name + std::string("|kmer_hash"));
  kmerge->add_dataset(counts_out, params.group_name + std::string("|count"));
    
  delete kmerge;

  leveldb::DB* db;
  leveldb::Options options;
  options.create_if_missing = true;

  leveldb::Status s = leveldb::DB::Open(options, params.db_filename, &db);
  REQUIRE(s.ok() == true);

  s = db->Get(leveldb::ReadOptions(), params.group_name + std::string("|size"), &value);
  REQUIRE(s.ok() == true);
  uint uncompressed_size = std::stoul(value);

  s = db->Get(leveldb::ReadOptions(), params.group_name + std::string("|kmer_hash"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(comp_hashes, ss_in);
  std::vector<uint> hashes_in = KMerge::uncompress(comp_hashes, uncompressed_size);
  comp_hashes.clear();
  std::vector<uint>().swap(comp_hashes);

  ss_in.str("");
  s = db->Get(leveldb::ReadOptions(), params.group_name + std::string("|count"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(comp_counts, ss_in);
  std::vector<uint> counts_in = KMerge::uncompress(comp_counts, uncompressed_size);
  comp_counts.clear();
  std::vector<uint>().swap(comp_counts);


  for (pos = 0; pos < hashes_in.size(); pos++) {
    REQUIRE(hashes_out[pos] == hashes_in[pos]);
    REQUIRE(counts_out[pos] == counts_in[pos]);
  }

  ss_delete << "rm " << params.db_filename << "/*";
  system(ss_delete.str().c_str());
  ss_delete.str("");
  ss_delete << "rm -r " << params.db_filename;
  system(ss_delete.str().c_str());
}


TEST_CASE("ThreadedParseKmerCountsAndCreateDBFromFastq", "[HashTest]") {
  param_struct params;
  std::string value;
  std::stringstream ss_in, ss_delete;
  std::vector<uint> comp_hashes, comp_counts;

  params.k_val_start = 3;
  params.k_val_end = 7;
  params.db_filename = "thread_fastq.db";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz";
  params.group_name = "sample";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  KMerge *kmerge = new KMerge(params.db_filename, "lookup3", ".");

  params.kmerge = kmerge;

  kmerge->build(params);

  delete kmerge;

  leveldb::DB* db;
  leveldb::Options options;
  options.create_if_missing = true;

  leveldb::Status s = leveldb::DB::Open(options, params.db_filename, &db);
  REQUIRE(s.ok() == true);

  s = db->Get(leveldb::ReadOptions(), params.group_name + std::string("|size"), &value);
  REQUIRE(s.ok() == true);
  uint uncompressed_size = std::stoul(value);

  REQUIRE(uncompressed_size == 8727);

  s = db->Get(leveldb::ReadOptions(), params.group_name + std::string("|kmer_hash"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(comp_hashes, ss_in);
  std::vector<uint> hashes_in = KMerge::uncompress(comp_hashes, uncompressed_size);
  comp_hashes.clear();
  std::vector<uint>().swap(comp_hashes);

  ss_in.str("");
  s = db->Get(leveldb::ReadOptions(), params.group_name + std::string("|count"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(comp_counts, ss_in);
  std::vector<uint> counts_in = KMerge::uncompress(comp_counts, uncompressed_size);
  comp_counts.clear();
  std::vector<uint>().swap(comp_counts);

  REQUIRE(hashes_in.size() == 8727);
  REQUIRE(counts_in.size() == 8727);

  delete db;

  ss_delete << "rm " << params.db_filename << "/*";
  system(ss_delete.str().c_str());
  ss_delete.str("");
  ss_delete << "rm -r " << params.db_filename;
  system(ss_delete.str().c_str());


}

TEST_CASE("ThreadedParseKmerCountsAndCreateDB", "[HashTest]") {
  std::vector<uint> hashes_in, counts_in, comp_hashes, comp_counts;
  uint kmer1_count, kmer2_count, kmer1_pos, kmer2_pos, pos, uncompressed_size;
  param_struct params1, params2, params3;
  std::map<std::string, std::string> taxonomy;
  std::string value, line;
  ifstream in_file;
  std::stringstream ss_in, ss_delete;

  const string kmer1("AAAAA");
  const string kmer2("GCGAT");
  leveldb::DB* db;
  leveldb::Options options;
  options.create_if_missing = true;


  params1.k_val_start = 5;
  params1.k_val_end = 5;
  params2.k_val_start = 5;
  params2.k_val_end = 5;
  params3.k_val_start = 5;
  params3.k_val_end = 5;

  params1.db_filename = "thread_example.db";
  params1.lock_filename = params1.db_filename + std::string(".lck");
  params2.db_filename = "thread_example.db";
  params2.lock_filename = params2.db_filename + std::string(".lck");
  params3.db_filename = "thread_example.db";
  params3.lock_filename = params3.db_filename + std::string(".lck");

  params1.seq_filename = "/home/darryl/Development/kmerge/tests/208831/208831.fasta.gz";
  params1.group_name = "208831";
  params1.num_threads = (params1.k_val_end - params1.k_val_start) / 2 + 1;
  params1.priority = 1;

  params2.seq_filename = "/home/darryl/Development/kmerge/tests/209328/209328.fasta.gz";
  params2.group_name = "209328";
  params2.num_threads = (params2.k_val_end - params2.k_val_start) / 2 + 1;
  params2.priority = 2;

  params3.seq_filename = "/home/darryl/Development/kmerge/tests/54095/54095.fasta.gz";
  params3.group_name = "54095";
  params3.num_threads = (params3.k_val_end - params3.k_val_start) / 2 + 1;
  params3.priority = 3;

  uint thread_count = 3;

  dlib::thread_pool tp(thread_count);

  KMerge* kmerge = new KMerge(params1.db_filename, "lookup3", ".");

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


  leveldb::Status s = leveldb::DB::Open(options, params1.db_filename, &db);
  REQUIRE(s.ok() == true);

  kmer1_count = 6150 /*AAAAA*/ + 6021 /*TTTTT*/;
  kmer2_count = 10775 /*GCGAT*/ + 10855 /*ATCGC*/;

  s = db->Get(leveldb::ReadOptions(), params1.group_name + std::string("|size"), &value);
  REQUIRE(s.ok() == true);
  uncompressed_size = std::stoul(value);


  s = db->Get(leveldb::ReadOptions(), params1.group_name + std::string("|kmer_hash"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(comp_hashes, ss_in);
  hashes_in = KMerge::uncompress(comp_hashes, uncompressed_size);
  comp_hashes.clear();
  std::vector<uint>().swap(comp_hashes);

  ss_in.str("");
  s = db->Get(leveldb::ReadOptions(), params1.group_name + std::string("|count"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(comp_counts, ss_in);
  counts_in = KMerge::uncompress(comp_counts, uncompressed_size);
  comp_counts.clear();
  std::vector<uint>().swap(comp_counts);

  for (pos = 0; pos < uncompressed_size; pos++) {
    if (hashes_in[pos] == KMerge::hash_kmer(kmer1, LOOKUP3)) {
      kmer1_pos = pos;
    }
    if (hashes_in[pos] == KMerge::hash_kmer(kmer2, LOOKUP3)) {
      kmer2_pos = pos;
    }
  }

  REQUIRE(hashes_in[kmer1_pos] == KMerge::hash_kmer(kmer1, LOOKUP3));
  REQUIRE(hashes_in[kmer2_pos] == KMerge::hash_kmer(kmer2, LOOKUP3));
  REQUIRE(counts_in[kmer1_pos] == kmer1_count);
  REQUIRE(counts_in[kmer2_pos] == kmer2_count);

  ss_in.str("");
  hashes_in.clear();
  std::vector<uint>().swap(hashes_in);
  counts_in.clear();
  std::vector<uint>().swap(counts_in);
 
 
  s = db->Get(leveldb::ReadOptions(), params1.group_name + std::string("|taxonomy"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(taxonomy, ss_in);


  in_file.open("/home/darryl/Development/kmerge/tests/208831/taxonomy.txt");

  while (std::getline(in_file, line)) {
    std::istringstream tokenizer(line);
    std::string classification, taxon;
    uint i = 0;
    while (!tokenizer.eof()) {
      std::string token;
      getline(tokenizer, token, '\t');
      if (i == 0) { // this is taxon                                                                                                                               
        taxon = token;
      } else { // this is classification                                                                                                                           
        classification = token;
      }
      i++;
    }
    REQUIRE(taxonomy[taxon] == classification);
  }
  taxonomy.clear();
  std::map<std::string, std::string>().swap(taxonomy);
  
  in_file.close();

  kmer1_count = 2147 /*AAAAA*/ + 1919 /*TTTTT*/;
  kmer2_count = 12082 /*GCGAT*/ + 12213 /*ATCGC*/;

  s = db->Get(leveldb::ReadOptions(), params2.group_name + std::string("|size"), &value);
  REQUIRE(s.ok() == true);
  uncompressed_size = std::stoul(value);


  s = db->Get(leveldb::ReadOptions(), params2.group_name + std::string("|kmer_hash"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(comp_hashes, ss_in);
  hashes_in = KMerge::uncompress(comp_hashes, uncompressed_size);
  comp_hashes.clear();
  std::vector<uint>().swap(comp_hashes);

  ss_in.str("");
  s = db->Get(leveldb::ReadOptions(), params2.group_name + std::string("|count"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(comp_counts, ss_in);
  counts_in = KMerge::uncompress(comp_counts, uncompressed_size);
  comp_counts.clear();
  std::vector<uint>().swap(comp_counts);

  for (pos = 0; pos < uncompressed_size; pos++) {
    if (hashes_in[pos] == KMerge::hash_kmer(kmer1, LOOKUP3)) {
      kmer1_pos = pos;
    }
    if (hashes_in[pos] == KMerge::hash_kmer(kmer2, LOOKUP3)) {
      kmer2_pos = pos;
    }
  }

  REQUIRE(hashes_in[kmer1_pos] == KMerge::hash_kmer(kmer1, LOOKUP3));
  REQUIRE(hashes_in[kmer2_pos] == KMerge::hash_kmer(kmer2, LOOKUP3));
  REQUIRE(counts_in[kmer1_pos] == kmer1_count);
  REQUIRE(counts_in[kmer2_pos] == kmer2_count);

  ss_in.str("");
  hashes_in.clear();
  std::vector<uint>().swap(hashes_in);
  counts_in.clear();
  std::vector<uint>().swap(counts_in);

  s = db->Get(leveldb::ReadOptions(), params2.group_name + std::string("|taxonomy"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(taxonomy, ss_in);

  in_file.open("/home/darryl/Development/kmerge/tests/209328/taxonomy.txt");

  while (std::getline(in_file, line)) {
    std::istringstream tokenizer(line);
    std::string classification, taxon;
    uint i = 0;
    while (!tokenizer.eof()) {
      std::string token;
      getline(tokenizer, token, '\t');
      if (i == 0) { // this is taxon                                                                                                                               
        taxon = token;
      } else { // this is classification                                                                                                                           
        classification = token;
      }
      i++;
    }
    REQUIRE(taxonomy[taxon] == classification);
  }
  taxonomy.clear();
  std::map<std::string, std::string>().swap(taxonomy);

  in_file.close();

  kmer1_count = 9896 /*AAAAA*/ + 9505 /*TTTTT*/;
  kmer2_count = 733 /*GCGAT*/ + 750 /*ATCGC*/;

  s = db->Get(leveldb::ReadOptions(), params3.group_name + std::string("|size"), &value);
  REQUIRE(s.ok() == true);
  uncompressed_size = std::stoul(value);


  s = db->Get(leveldb::ReadOptions(), params3.group_name + std::string("|kmer_hash"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(comp_hashes, ss_in);
  hashes_in = KMerge::uncompress(comp_hashes, uncompressed_size);
  comp_hashes.clear();
  std::vector<uint>().swap(comp_hashes);

  ss_in.str("");
  s = db->Get(leveldb::ReadOptions(), params3.group_name + std::string("|count"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(comp_counts, ss_in);
  counts_in = KMerge::uncompress(comp_counts, uncompressed_size);
  comp_counts.clear();
  std::vector<uint>().swap(comp_counts);

  for (pos = 0; pos < uncompressed_size; pos++) {
    if (hashes_in[pos] == KMerge::hash_kmer(kmer1, LOOKUP3)) {
      kmer1_pos = pos;
    }
    if (hashes_in[pos] == KMerge::hash_kmer(kmer2, LOOKUP3)) {
      kmer2_pos = pos;
    }
  }

  REQUIRE(hashes_in[kmer1_pos] == KMerge::hash_kmer(kmer1, LOOKUP3));
  REQUIRE(hashes_in[kmer2_pos] == KMerge::hash_kmer(kmer2, LOOKUP3));
  REQUIRE(counts_in[kmer1_pos] == kmer1_count);
  REQUIRE(counts_in[kmer2_pos] == kmer2_count);

  ss_in.str("");
  hashes_in.clear();
  std::vector<uint>().swap(hashes_in);
  counts_in.clear();
  std::vector<uint>().swap(counts_in);

  s = db->Get(leveldb::ReadOptions(), params3.group_name + std::string("|taxonomy"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(taxonomy, ss_in);


  in_file.open("/home/darryl/Development/kmerge/tests/54095/taxonomy.txt");

  while (std::getline(in_file, line)) {
    std::istringstream tokenizer(line);
    std::string classification, taxon;
    uint i = 0;
    while (!tokenizer.eof()) {
      std::string token;
      getline(tokenizer, token, '\t');
      if (i == 0) { // this is taxon                                                                                                                               
        taxon = token;
      } else { // this is classification                                                                                                                           
        classification = token;
      }
      i++;
    }
    REQUIRE(taxonomy[taxon] == classification);
  }
  taxonomy.clear();
  std::map<std::string, std::string>().swap(taxonomy);

  in_file.close();

  delete db;

  ss_delete << "rm " << params1.db_filename << "/*";
  system(ss_delete.str().c_str());
  ss_delete.str("");
  ss_delete << "rm -r " << params1.db_filename;
  system(ss_delete.str().c_str());

}

TEST_CASE("TestHashingFunctions", "[HashTest]") {
  ulib::chain_hash_map<uint, uint> hashed_counts(100000000);
  ulib::chain_hash_map<uint, uint>::const_iterator m_iter;
  param_struct params;

  params.k_val_start = 5;
  params.k_val_end = 5;
  params.seq_filename = "/home/darryl/Development/kmerge/tests/genome.test.fasta.gz";
  params.group_name = "org1";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  params.kmerge = new KMerge(params.db_filename, "lookup3", ".");
  bool success = params.kmerge->count_hashed_kmers(params, hashed_counts, true);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter.key() >= 0);
    REQUIRE(m_iter.value() <= (uint) MAX_UINT_VAL);
  }
  hashed_counts.clear();
  delete params.kmerge;

  params.kmerge = new KMerge(params.db_filename, "spooky", ".");
  success = params.kmerge->count_hashed_kmers(params, hashed_counts, true);
  REQUIRE(success == true);

  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter.key() >= 0);
    REQUIRE(m_iter.value() <= (uint) MAX_UINT_VAL);
  }
  hashed_counts.clear();
  delete params.kmerge;

  params.kmerge = new KMerge(params.db_filename, "city", ".");
  success = params.kmerge->count_hashed_kmers(params, hashed_counts, true);
  REQUIRE(success == true);


  for(m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter.key() >= 0);
    REQUIRE(m_iter.value() <= (uint) MAX_UINT_VAL);
  }
  hashed_counts.clear();
  delete params.kmerge;

}


TEST_CASE("AddTaxonomyInfoToDB", "[LevelDBTest]") {
  std::string db_filename("taxonomy.db"), group("54095"); 
  std::string line, value;
  std::stringstream ss_in, ss_delete;
  std::map<std::string, std::string> taxonomy;

  

  KMerge *kmerge = new KMerge(db_filename, "lookup3", ".");

  kmerge->add_taxonomy(group);

  delete kmerge;

  leveldb::DB* db;
  leveldb::Options options;
  options.create_if_missing = true;

  leveldb::Status s = leveldb::DB::Open(options, db_filename, &db);
  REQUIRE(s.ok() == true);

   
  s = db->Get(leveldb::ReadOptions(), group + std::string("|taxonomy"), &value);
  REQUIRE(s.ok() == true);
  ss_in << value;
  dlib::deserialize(taxonomy, ss_in);


  ifstream in_file("/home/darryl/Development/kmerge/tests/54095/taxonomy.txt");

  while (std::getline(in_file, line)) {
    std::istringstream tokenizer(line);
    std::string classification, taxon;
    uint i = 0;
    while (!tokenizer.eof()) {
      std::string token;
      getline(tokenizer, token, '\t');
      if (i == 0) { // this is taxon                                                                                                                               
        taxon = token;
      } else { // this is classification                                                                                                                           
        classification = token;
      }
      i++;
    }
    REQUIRE(taxonomy[taxon] == classification);
  }
  taxonomy.clear();
  std::map<std::string, std::string>().swap(taxonomy);

  delete db;

  
  ss_delete << "rm " << db_filename << "/*";
  system(ss_delete.str().c_str());
  ss_delete.str("");
  ss_delete << "rm -r " << db_filename;
  system(ss_delete.str().c_str());


}

