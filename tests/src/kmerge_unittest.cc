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
#include <sys/stat.h>



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
  btree::btree_map<std::string, uint> kmer_counts, jf_counts;
  kseq_t *seq;
  gzFile fp;
  std::vector<std::string> kmers;

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
      kmers.push_back(kmer);
    }
  }
  kseq_destroy(seq);
  gzclose(fp);
  
  for (auto it = kmers.begin(); it != kmers.end(); it++) {
    kmer_counts[*it]++;
  }

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
  std::vector<std::string> kmers;

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
      kmers.push_back(kmer);
    }
  }
  kseq_destroy(seq);
  gzclose(fp);

  for (auto it = kmers.begin(); it != kmers.end(); it++) {
    kmer_counts[*it]++;
  }

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
  btree::btree_map<uint,uint> hashed_counts;
  param_struct params;
  std::vector<std::string> files;
  kseq_t *seq;
  gzFile fp;
  int l;
  std::mutex mtx;


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

  KMerge *kmerge = new KMerge(params.db_filename.c_str(), "lookup3", ".", MIN_MEM);

  params.kmerge = kmerge;
 
  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s);
    for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
      bool success = params.kmerge->hash_seq(seq_str, k, hashed_counts, mtx);
      REQUIRE(success == true);
    }
  }
  kseq_destroy(seq);
  gzclose(fp);


  std::string seq_str2, last = "";
  uint count;
  for (std::string filename : files ) {
    ifstream in_file (filename.c_str());
    if (in_file.is_open()) {
      while ( !in_file.eof() ) {
	in_file >> seq_str2 >> count;
	if (seq_str2 != last) {
	  uint hash = kmerge->hash_kmer(seq_str2);
	  jf_hashed_counts[hash] = jf_hashed_counts[hash] +  count;
	}
	last = seq_str2;
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
  btree::btree_map<uint, uint> hashed_counts, jf_hashed_counts;
  param_struct params;
  std::vector<std::string> files;
  kseq_t *seq;
  gzFile fp;
  int l;
  std::mutex mtx;


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


  KMerge *kmerge = new KMerge(params.db_filename.c_str(), "lookup3", ".", MIN_MEM);

  params.kmerge = kmerge;

  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s);
    for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
      bool success = params.kmerge->hash_seq(seq_str, k, hashed_counts, mtx);
      REQUIRE(success == true);
    }
  }
  kseq_destroy(seq);
  gzclose(fp);



  std::string seq_str2, last = "";
  uint count;
  for (std::string filename : files ) {
    ifstream in_file (filename.c_str());
    if (in_file.is_open()) {
      while ( !in_file.eof() ) {
	in_file >> seq_str2 >> count;
	if (seq_str2 != last) {
	  uint hash = kmerge->hash_kmer(seq_str2);
	  jf_hashed_counts[hash] = jf_hashed_counts[hash] +  count;
	}
	last = seq_str2;
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
  btree::btree_map<uint, uint> hashed_counts, hashed_counts2;
  int l;
  param_struct params;
  std::vector<std::tuple<uint, uint, uint, uint> > coords;
  uint pieces, piece_length;
  std::mutex mtx;
  kseq_t *seq;
  gzFile fp;
  std::vector<std::string> seqs;
  uint seq_id;
  
  params.k_val_start = 3;
  params.k_val_end = 11;
  params.db_filename = "dummy.db";
  params.seq_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/208831/208831.fasta.gz";
  params.group_name = "208831";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  KMerge *kmerge = new KMerge(params.db_filename, "lookup3", ".", MIN_MEM);
  params.kmerge = kmerge;

 
  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  seq_id = 0;
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s);
    seqs.push_back(seq_str);
    uint str_len = seq_str.size();
    for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
      if (str_len < k) {
	continue;
      } else {
	pieces = params.num_threads;
	piece_length = str_len / pieces;
      }
      uint pos = 0;
      while (pos < str_len) {
	uint start = pos, end = pos + piece_length + k - 1;
	if (end > str_len) end = str_len;
	coords.push_back(std::make_tuple(seq_id, k, start, end));
	pos += piece_length;
	if (end == str_len) break; // all sequence has been accounted for
      }
    } 
    seq_id++;
  }
  KMerge::HashSeq func(params, seqs, hashed_counts, coords, mtx, false);
  dlib::parallel_for(params.num_threads, 0, coords.size(), func);
  coords.clear();
  kseq_destroy(seq);
  gzclose(fp);
  

 
  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s);
    for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
      bool success = params.kmerge->hash_seq(seq_str, k, hashed_counts2, mtx);
      REQUIRE(success == true);
    }
  }
  kseq_destroy(seq);
  gzclose(fp);

  REQUIRE(hashed_counts.size() == hashed_counts2.size());

  for (auto it = hashed_counts.begin(); it != hashed_counts.end(); it++) {
    REQUIRE(hashed_counts[it->first] == hashed_counts2[it->first]);
  }

  delete kmerge;
}

TEST_CASE("CountHashedKmersInParallelFastq", "[HashTest]") {
  btree::btree_map<uint, uint> hashed_counts, hashed_counts2;
  int l;
  param_struct params;
  std::vector<std::tuple<uint, uint, uint, uint> > coords;
  uint pieces, piece_length;
  std::mutex mtx;
  kseq_t *seq;
  gzFile fp;
  std::vector<std::string> seqs;
  uint seq_id;
  
  params.k_val_start = 3;
  params.k_val_end = 11;
  params.db_filename = "dummy.db";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz";
  params.group_name = "sample";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  KMerge *kmerge = new KMerge(params.db_filename, "lookup3", ".", MIN_MEM);
  params.kmerge = kmerge;


  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  seq_id = 0;
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s);
    seqs.push_back(seq_str);
    uint str_len = seq_str.size();
    for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
      if (str_len < k) {
	continue;
      } else {
	pieces = params.num_threads;
	piece_length = str_len / pieces;
      }
      uint pos = 0;
      while (pos < str_len) {
	uint start = pos, end = pos + piece_length + k - 1;
	if (end > str_len) end = str_len;
	coords.push_back(std::make_tuple(seq_id, k, start, end));
	pos += piece_length;
	if (end == str_len) break; // all sequence has been accounted for
      }
    } 
    seq_id++;
  }
  KMerge::HashSeq func(params, seqs, hashed_counts, coords, mtx, false);
  dlib::parallel_for(params.num_threads, 0, coords.size(), func);
  coords.clear();
  kseq_destroy(seq);
  gzclose(fp);
  

 
  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s);
    for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
      bool success = params.kmerge->hash_seq(seq_str, k, hashed_counts2, mtx);
      REQUIRE(success == true);
    }
  }
  kseq_destroy(seq);
  gzclose(fp);

  REQUIRE(hashed_counts.size() == hashed_counts2.size());

  for (auto it = hashed_counts.begin(); it != hashed_counts.end(); it++) {
    REQUIRE(hashed_counts[it->first] == hashed_counts2[it->first]);
  }

  delete kmerge;
}


TEST_CASE("ParseKmerCountsAndCreateDB", "[HashTest]") {
  std::vector<uint>comp_counts;
  ulib::chain_hash_map<uint, uint> hashed_counts(100000000);
  const string kmer1("AAAAA");
  const string kmer2("TTTTT");
  const uint kmer1_count = 84;
  const uint kmer2_count = 13;
  uint kmer1_pos = 0;
  uint kmer2_pos = 0;
  uint pos = 0;
  param_struct params;
  dlib::thread_pool tp(1);
 
  params.k_val_start = 5;
  params.k_val_end = 13;
  params.db_filename = "/home/darryl/Development/kmerge/tests/dump_example.db";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/208831/208831.fasta.gz";
  params.group_name = "208831";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
  params.priority = 1;
  params.lock_filename = params.db_filename + std::string(".lck");
  params.is_ref = true;

  KMerge* kmerge = new KMerge(params.db_filename, "lookup3", ".", 0.5);
  params.kmerge = kmerge;
  
  KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
  tp.add_task(*task, &KMerge::BuilderTask::execute);

  tp.wait_for_all_tasks();
 
  delete kmerge;
 
  int rc = unqlite_open(&(params.db), params.db_filename.c_str(), UNQLITE_OPEN_CREATE);
  REQUIRE(rc == UNQLITE_OK);


  unqlite_int64 n_bytes;
  n_bytes = sizeof(uint);
  uint size;
  rc = unqlite_kv_fetch(params.db,(params.group_name + std::string("|kmer_hash|0|size")).c_str(),-1, &size, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);

  uint *value = new uint[size];
  n_bytes = size*sizeof(uint);
  rc = unqlite_kv_fetch(params.db,(params.group_name + std::string("|kmer_hash|0")).c_str(),-1,value, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);
  std::vector<uint> hashes_in;
  hashes_in.assign((uint*) &value[0], (uint*) &value[0] + size);
  delete [] value;
 
  n_bytes = sizeof(uint);
  uint compressed_size;
  rc = unqlite_kv_fetch(params.db,(params.group_name + std::string("|count|0|size")).c_str(),-1, &compressed_size, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);
  value = new uint[compressed_size];
  n_bytes = sizeof(uint)*compressed_size;
  rc = unqlite_kv_fetch(params.db,(params.group_name + std::string("|count|0")).c_str(),-1,value, &n_bytes);
  comp_counts.assign((uint*) &value[0], (uint*) &value[0] + compressed_size);
  std::vector<uint> counts_in = KMerge::uncompress(comp_counts, size);
  comp_counts.clear();
  std::vector<uint>().swap(comp_counts);
  delete [] value;

  unqlite_close(params.db);

  std::partial_sum(hashes_in.begin(), hashes_in.end(), hashes_in.begin());

  REQUIRE(hashes_in.size() == 6091121);
  REQUIRE(counts_in.size() == 6091121);


  if( remove( params.db_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

}


TEST_CASE("ThreadedParseKmerCountsAndCreateDBFromFastq", "[HashTest]") {
  param_struct params;
  std::string value;
  std::stringstream ss_delete;
  std::vector<uint> comp_hashes;
  std::vector<uint> comp_counts;
  dlib::thread_pool tp(5);
  std::mutex mtx;

  params.k_val_start = 3;
  params.k_val_end = 7;
  params.db_filename = "thread_fq.db";
  params.lock_filename = params.db_filename + std::string(".lck");
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz";
  params.group_name = "sample";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
  params.is_ref = false;

  KMerge *kmerge = new KMerge(params.db_filename, "lookup3", ".", MIN_MEM);

  params.kmerge = kmerge;

  KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
  tp.add_task(*task, &KMerge::BuilderTask::execute);
  tp.wait_for_all_tasks();



  delete kmerge;

  unqlite *db;
  int rc = unqlite_open(&db, params.db_filename.c_str(), UNQLITE_OPEN_CREATE);
  REQUIRE(rc == UNQLITE_OK);

  unqlite_int64 n_bytes;
  n_bytes = sizeof(uint);
  uint num_parts;
  rc = unqlite_kv_fetch(db, (params.group_name + std::string("|parts")).c_str(), -1, &num_parts, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);

  std::vector<uint> hashes_in;
  std::vector<uint> counts_in;

  for (uint i=0; i < num_parts; i++) {
    n_bytes = sizeof(uint);
    uint size;
    rc = unqlite_kv_fetch(db,(params.group_name + std::string("|kmer_hash|") + std::to_string(i) + std::string("|size")).c_str(),-1, &size, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);

    REQUIRE(size == 8727);

    uint *value = new uint[size];
    n_bytes = size*sizeof(uint);
    rc = unqlite_kv_fetch(db,(params.group_name + std::string("|kmer_hash|") + std::to_string(i)).c_str(),-1,value, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);
    std::vector<uint> part_hashes;
    part_hashes.assign((uint*) &value[0], (uint*) &value[0] + size);
    delete [] value;

    n_bytes = sizeof(uint);
    uint compressed_size;
    rc = unqlite_kv_fetch(db,(params.group_name + std::string("|count|") + std::to_string(i) + std::string("|size")).c_str(),-1, &compressed_size, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);
    value = new uint[compressed_size];
    n_bytes = sizeof(uint)*compressed_size;
    rc = unqlite_kv_fetch(db,(params.group_name + std::string("|count|") + std::to_string(i)).c_str(),-1,value, &n_bytes);
    comp_counts.assign((uint*) &value[0], (uint*) &value[0] + compressed_size);
    std::vector<uint> part_counts = KMerge::uncompress(comp_counts, size);
    comp_counts.clear();
    std::vector<uint>().swap(comp_counts);
    delete [] value;

    
    hashes_in.insert(hashes_in.end(), part_hashes.begin(), part_hashes.end());
    counts_in.insert(counts_in.end(), part_counts.begin(), part_counts.end());
  }
 
  std::partial_sum(hashes_in.begin(), hashes_in.end(), hashes_in.begin());
  
  REQUIRE(hashes_in.size() == 8727);
  REQUIRE(counts_in.size() == 8727);
  
  //ensure that hashes are in sorted order
  uint last = 0;
  for (std::vector<uint>::const_iterator v_iter = hashes_in.begin(); v_iter != hashes_in.end(); v_iter++) {
    if (v_iter != hashes_in.begin()) {
      REQUIRE(*v_iter > last);
    }
    last = *v_iter;
  }

  unqlite_close(db);

  if( remove( params.db_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

}



TEST_CASE("ThreadedParseKmerCountsAndCreateDB", "[HashTest]") {
  std::vector<uint> hashes_in, comp_hashes;
  std::vector<uint> counts_in, comp_counts;
  uint kmer1_count, kmer2_count, kmer1_pos, kmer2_pos, pos, total_uncompressed_size;
  param_struct params1, params2, params3;
  std::map<std::string, std::string> taxonomy;
  std::string line;
  char* value;
  ifstream in_file;
  std::stringstream ss_in, ss_delete;
  unqlite_int64 n_bytes;
  uint num_parts;

  const string kmer1("AAAAA");
  const string kmer2("GCGAT");


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
  params1.is_ref = true;

  params2.seq_filename = "/home/darryl/Development/kmerge/tests/209328/209328.fasta.gz";
  params2.group_name = "209328";
  params2.num_threads = (params2.k_val_end - params2.k_val_start) / 2 + 1;
  params2.priority = 2;
  params2.is_ref = true;

  params3.seq_filename = "/home/darryl/Development/kmerge/tests/54095/54095.fasta.gz";
  params3.group_name = "54095";
  params3.num_threads = (params3.k_val_end - params3.k_val_start) / 2 + 1;
  params3.priority = 3;
  params3.is_ref = true;

  uint thread_count = 3;

  dlib::thread_pool tp(thread_count);



  KMerge* kmerge = new KMerge(params1.db_filename, "lookup3", ".", 0.5);

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

  unqlite *db;
  int rc = unqlite_open(&db, params1.db_filename.c_str(), UNQLITE_OPEN_CREATE);
  REQUIRE(rc == UNQLITE_OK);

  kmer1_count = 6150 /*AAAAA*/ + 6021 /*TTTTT*/;
  kmer2_count = 10775 /*GCGAT*/ + 10855 /*ATCGC*/;

  n_bytes = sizeof(uint);
  rc = unqlite_kv_fetch(db, (params1.group_name + std::string("|parts")).c_str(), -1, &num_parts, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);
  for (uint i=0; i < num_parts; i++) {
    n_bytes = sizeof(uint);
    uint size;
    rc = unqlite_kv_fetch(db,(params1.group_name + std::string("|kmer_hash|") + std::to_string(i) + std::string("|size")).c_str(),-1, &size, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);

    uint *value = new uint[size];
    n_bytes = size*sizeof(uint);
    rc = unqlite_kv_fetch(db,(params1.group_name + std::string("|kmer_hash|") + std::to_string(i)).c_str(),-1,value, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);
    std::vector<uint> part_hashes;
    part_hashes.assign((uint*) &value[0], (uint*) &value[0] + size);
    delete [] value;

    n_bytes = sizeof(uint);
    uint compressed_size;
    rc = unqlite_kv_fetch(db,(params1.group_name + std::string("|count|") + std::to_string(i) + std::string("|size")).c_str(),-1, &compressed_size, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);
    value = new uint[compressed_size];
    n_bytes = sizeof(uint)*compressed_size;
    rc = unqlite_kv_fetch(db,(params1.group_name + std::string("|count|") + std::to_string(i)).c_str(),-1,value, &n_bytes);
    comp_counts.assign((uint*) &value[0], (uint*) &value[0] + compressed_size);
    std::vector<uint> part_counts = KMerge::uncompress(comp_counts, size);
    comp_counts.clear();
    std::vector<uint>().swap(comp_counts);
    delete [] value;

    
    hashes_in.insert(hashes_in.end(), part_hashes.begin(), part_hashes.end());
    counts_in.insert(counts_in.end(), part_counts.begin(), part_counts.end());
    
  }

  std::partial_sum(hashes_in.begin(), hashes_in.end(), hashes_in.begin());

  REQUIRE(hashes_in.size() == 512);
  REQUIRE(counts_in.size() == 512);

  //ensure that hashes are in sorted order
  uint last = 0;
  for (std::vector<uint>::const_iterator v_iter = hashes_in.begin(); v_iter != hashes_in.end(); v_iter++) {
    if (v_iter != hashes_in.begin()) {
      REQUIRE(*v_iter > last);
    }
    last = *v_iter;
  }

  for (pos = 0; pos < 512 /*3^5*/; pos++) {
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

  hashes_in.clear();
  std::vector<uint>().swap(hashes_in);
  counts_in.clear();
  std::vector<uint>().swap(counts_in);
 
  n_bytes = sizeof(uint);
  uint size;
  rc = unqlite_kv_fetch(db,(params1.group_name + std::string("|taxonomy") + std::string("|size")).c_str(),-1, &size, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);


  n_bytes = size;
  char* tax_value = new char[size];
  rc = unqlite_kv_fetch(db,(params1.group_name + std::string("|taxonomy")).c_str(),-1, tax_value, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);
  ss_in << tax_value;
  dlib::deserialize(taxonomy, ss_in);
  delete [] tax_value;

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

  n_bytes = sizeof(uint);
  rc = unqlite_kv_fetch(db, (params2.group_name + std::string("|parts")).c_str(), -1, &num_parts, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);

  for (uint i=0; i < num_parts; i++) {
    n_bytes = sizeof(uint);
    uint size;
    rc = unqlite_kv_fetch(db,(params2.group_name + std::string("|kmer_hash|") + std::to_string(i) + std::string("|size")).c_str(),-1, &size, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);

    uint *value = new uint[size];
    n_bytes = size*sizeof(uint);
    rc = unqlite_kv_fetch(db,(params2.group_name + std::string("|kmer_hash|") + std::to_string(i)).c_str(),-1,value, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);
    std::vector<uint> part_hashes;
    part_hashes.assign((uint*) &value[0], (uint*) &value[0] + size);
    delete [] value;

    n_bytes = sizeof(uint);
    uint compressed_size;
    rc = unqlite_kv_fetch(db,(params2.group_name + std::string("|count|") + std::to_string(i) + std::string("|size")).c_str(),-1, &compressed_size, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);
    value = new uint[compressed_size];
    n_bytes = sizeof(uint)*compressed_size;
    rc = unqlite_kv_fetch(db,(params2.group_name + std::string("|count|") + std::to_string(i)).c_str(),-1,value, &n_bytes);
    comp_counts.assign((uint*) &value[0], (uint*) &value[0] + compressed_size);
    std::vector<uint> part_counts = KMerge::uncompress(comp_counts, size);
    comp_counts.clear();
    std::vector<uint>().swap(comp_counts);
    delete [] value;

    
    hashes_in.insert(hashes_in.end(), part_hashes.begin(), part_hashes.end());
    counts_in.insert(counts_in.end(), part_counts.begin(), part_counts.end());
    
  }

  std::partial_sum(hashes_in.begin(), hashes_in.end(), hashes_in.begin());

  //ensure that hashes are in sorted order
  for (std::vector<uint>::const_iterator v_iter = hashes_in.begin(); v_iter != hashes_in.end(); v_iter++) {
    if (v_iter != hashes_in.begin()) {
      REQUIRE(*v_iter > last);
    }
    last = *v_iter;
  }

  for (pos = 0; pos < 512 /*3^5*/; pos++) {
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

  hashes_in.clear();
  std::vector<uint>().swap(hashes_in);
  counts_in.clear();
  std::vector<uint>().swap(counts_in);

  n_bytes = sizeof(uint);
  rc = unqlite_kv_fetch(db,(params2.group_name + std::string("|taxonomy") + std::string("|size")).c_str(),-1, &size, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);

  ss_in.str("");
  n_bytes = size;
  tax_value = new char[size];
  rc = unqlite_kv_fetch(db,(params2.group_name + std::string("|taxonomy")).c_str(),-1, tax_value, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);
  ss_in << tax_value;
  dlib::deserialize(taxonomy, ss_in);
  delete [] tax_value;


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

  n_bytes = sizeof(uint);
  rc = unqlite_kv_fetch(db, (params3.group_name + std::string("|parts")).c_str(), -1, &num_parts, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);

  for (uint i=0; i < num_parts; i++) {
    n_bytes = sizeof(uint);
    uint size;
    rc = unqlite_kv_fetch(db,(params3.group_name + std::string("|kmer_hash|") + std::to_string(i) + std::string("|size")).c_str(),-1, &size, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);

    uint *value = new uint[size];
    n_bytes = size*sizeof(uint);
    rc = unqlite_kv_fetch(db,(params3.group_name + std::string("|kmer_hash|") + std::to_string(i)).c_str(),-1,value, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);
    std::vector<uint> part_hashes;
    part_hashes.assign((uint*) &value[0], (uint*) &value[0] + size);
    delete [] value;

    n_bytes = sizeof(uint);
    uint compressed_size;
    rc = unqlite_kv_fetch(db,(params3.group_name + std::string("|count|") + std::to_string(i) + std::string("|size")).c_str(),-1, &compressed_size, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);
    value = new uint[compressed_size];
    n_bytes = sizeof(uint)*compressed_size;
    rc = unqlite_kv_fetch(db,(params3.group_name + std::string("|count|") + std::to_string(i)).c_str(),-1,value, &n_bytes);
    comp_counts.assign((uint*) &value[0], (uint*) &value[0] + compressed_size);
    std::vector<uint> part_counts = KMerge::uncompress(comp_counts, size);
    comp_counts.clear();
    std::vector<uint>().swap(comp_counts);
    delete [] value;

    
    hashes_in.insert(hashes_in.end(), part_hashes.begin(), part_hashes.end());
    counts_in.insert(counts_in.end(), part_counts.begin(), part_counts.end());
  }

  std::partial_sum(hashes_in.begin(), hashes_in.end(), hashes_in.begin());

  //ensure that hashes are in sorted order
  for (std::vector<uint>::const_iterator v_iter = hashes_in.begin(); v_iter != hashes_in.end(); v_iter++) {
    if (v_iter != hashes_in.begin()) {
      REQUIRE(*v_iter > last);
    }
    last = *v_iter;
  }


  for (pos = 0; pos < 512 /*3^5*/; pos++) {
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

  hashes_in.clear();
  std::vector<uint>().swap(hashes_in);
  counts_in.clear();
  std::vector<uint>().swap(counts_in);

  n_bytes = sizeof(uint);
  rc = unqlite_kv_fetch(db,(params3.group_name + std::string("|taxonomy") + std::string("|size")).c_str(),-1, &size, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);

  ss_in.str("");
  n_bytes = size;
  tax_value = new char[size];
  rc = unqlite_kv_fetch(db,(params3.group_name + std::string("|taxonomy")).c_str(),-1, tax_value, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);
  ss_in << tax_value;
  dlib::deserialize(taxonomy, ss_in);
  delete [] tax_value;



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

  unqlite_close(db);

  if( remove( params1.db_filename.c_str() ) != 0 )
    perror( "Error deleting file" );


}

TEST_CASE("TestHashingFunctions", "[HashTest]") {
  btree::btree_map<uint, uint> hashed_counts;
  param_struct params;

  params.k_val_start = 5;
  params.k_val_end = 5;
  params.seq_filename = "/home/darryl/Development/kmerge/tests/genome.test.fasta.gz";
  params.group_name = "org1";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
  params.is_ref = true;


  params.kmerge = new KMerge(params.db_filename, "lookup3", ".", MIN_MEM);
  bool success = params.kmerge->count_hashed_kmers(params, hashed_counts, true);
  REQUIRE(success == true);

  for(auto m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->second <= (uint) MAX_UINT_VAL);
  }
  hashed_counts.clear();
  delete params.kmerge;

  params.kmerge = new KMerge(params.db_filename, "spooky", ".", MIN_MEM);
  success = params.kmerge->count_hashed_kmers(params, hashed_counts, true);
  REQUIRE(success == true);

  for(auto m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->second <= (uint) MAX_UINT_VAL);
  }
  hashed_counts.clear();
  delete params.kmerge;

  params.kmerge = new KMerge(params.db_filename, "city", ".", MIN_MEM);
  success = params.kmerge->count_hashed_kmers(params, hashed_counts, true);
  REQUIRE(success == true);


  for(auto m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->second <= (uint) MAX_UINT_VAL);
  }
  hashed_counts.clear();
  delete params.kmerge;

}


TEST_CASE("AddTaxonomyInfoToDB", "[LevelDBTest]") {
  std::string db_filename("taxonomy.db"), group("54095"); 
  std::string line;
  char* value;
  std::stringstream ss_in;
  std::map<std::string, std::string> taxonomy;

  unqlite *db;

  int rc = unqlite_open(&db, db_filename.c_str(), UNQLITE_OPEN_CREATE);
  REQUIRE(rc == UNQLITE_OK);

  
  KMerge *kmerge = new KMerge(db_filename, "lookup3", ".", MIN_MEM);

  kmerge->add_taxonomy(group, db);

  delete kmerge;

  unqlite_int64 n_bytes = sizeof(uint);
  uint size;
  rc = unqlite_kv_fetch(db,(group + std::string("|taxonomy") + std::string("|size")).c_str(),-1, &size, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);

  n_bytes = size;
  value = new char[size];
  rc = unqlite_kv_fetch(db,(group + std::string("|taxonomy")).c_str(),-1, value, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);
  ss_in << value;
  dlib::deserialize(taxonomy, ss_in);
  delete [] value;


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

  unqlite_close(db);

  if( remove( db_filename.c_str() ) != 0 )
    perror( "Error deleting file" );


}

