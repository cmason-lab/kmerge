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
#include <chrono>
#include <random>
#include <google/sparsetable>

TEST_CASE("CompressHashesTest", "CompressionTest") {
  size_t N = 10 * 1000;
  std::vector<uint> mydata(N);
  for(uint32_t i = 0; i < N;i += 150) mydata[i] = i;

  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > compressed_output = KMerge::compress(mydata);

  std::cout<<std::setprecision(3);
  std::cout<<"You are using " << 32.0 * static_cast<double>(compressed_output.size()) /
    static_cast<double>(mydata.size()) <<" bits per integer. "<<std::endl;
  

  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > mydataback = KMerge::uncompress(compressed_output, N);

  for (uint i = 0; i <= mydata.size(); i++) {
    REQUIRE(mydata[i] == mydataback[i]);
  }

}


TEST_CASE("UnQLiteBasicsTest", "HashStorageTest") {
  std::string key1("test"), value;
  uint value_out = 1, value_in;
  unqlite *db;
  int rc;
  unqlite_int64 n_bytes = sizeof(uint);

  rc = unqlite_open(&db,"test.db", UNQLITE_OPEN_CREATE);
  REQUIRE(rc == UNQLITE_OK);

  rc = unqlite_kv_store(db,key1.c_str(),-1,&value_out,sizeof(uint));
  REQUIRE(rc == UNQLITE_OK);

  rc = unqlite_kv_fetch(db,key1.c_str(),-1,&value_in, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);

  REQUIRE(value_in == 1);
  unqlite_close(db);

  if( remove( "test.db" ) != 0 )
    perror( "Error deleting file" );

}

TEST_CASE("UnQLiteAndFastPForTest", "HashStorageTest") {
  std::string key1("test");
  std::stringstream ss_in, ss_out;
  unqlite *db;
  int rc;
  unqlite_int64 n_bytes;
  uint * value;

  std::cout << "Generating original data" << std::endl;

  uint N = 4000000050;
  size_t step = 1500;
  //std::vector<uint> test;
  //REQUIRE(test.max_size() > N);
  std::vector<uint> mydata(N);
  for(uint32_t i = 0; i < N;i += step) mydata[i] = i;

  std::cout << "Finished generating original data" << std::endl;


  rc = unqlite_open(&db,"test.db", UNQLITE_OPEN_CREATE);
  REQUIRE(rc == UNQLITE_OK);

  std::cout << "Compressing data" << std::endl;
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > compressed = KMerge::compress(mydata);
  std::cout << "Finished compressing data" << std::endl;    

  std::cout << "Writing data" << std::endl;
  uint comp_size = compressed.size();
  rc = unqlite_kv_store(db,key1.c_str(),-1,&compressed[0],sizeof(uint)*compressed.size());
  REQUIRE(rc == UNQLITE_OK);
  std::cout << "Finished writing data" << std::endl; 
  compressed.clear();
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >().swap(compressed);

  std::cout << "Reading data" << std::endl;
  n_bytes = sizeof(uint)*comp_size;
  value = new uint [comp_size];
  rc = unqlite_kv_fetch(db,key1.c_str(),-1, value, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);

  std::cout << "Finished reading data" << std::endl;
  
  std::cout << "Copying data" << std::endl;
  compressed.assign((uint*) &value[0], (uint*) &value[0] + comp_size);
  std::cout << "Finished copying data" << std::endl;
  delete [] value;
  std::cout << "Uncompressing data" << std::endl;
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > mydataback = KMerge::uncompress(compressed, N);
  std::cout << "Finished uncompressing data" << std::endl;

  std::cout << "Checking data" << std::endl;
  REQUIRE(mydata.size() == mydataback.size());

  for (uint i = 0; i <= mydataback.size(); i+=step) {
    REQUIRE(mydataback[i] == i);
  }
  std::cout << "Finished checking data" << std::endl;

  unqlite_close(db);

  if( remove( "test.db" ) != 0 )
    perror( "Error deleting file" );

}


TEST_CASE("UnQLiteAndFastPForWithLargeDenseHashCountTest", "HashStorageTest") {
  std::string key1("test");
  std::stringstream ss_in, ss_out;
  unqlite *db;
  int rc;
  std::vector<uint> comp_size;
  unqlite_int64 n_bytes;
  uint * value;


  std::cout << "Generating original data" << std::endl;

  size_t N = 1352554547; // hash count causing failure in live runs
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > mydata(N);
  for(uint32_t i = 0; i < N; i++) mydata[i] = i;
  uint part_size = 500000000;
  std::cout << "Finished generating original data" << std::endl;

  rc = unqlite_open(&db,"test.db", UNQLITE_OPEN_CREATE);
  REQUIRE(rc == UNQLITE_OK);

  for (uint i=0; i < ceil((double)N/part_size); i++) {
    uint start_offset = i*part_size;
    uint end_offset = (i+1)*part_size;
    if (N - i*part_size < part_size) end_offset = N;
    std::vector<uint> sub(mydata.begin() + start_offset, mydata.begin() + end_offset);
    std::cout << "Compressing data" << std::endl;
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > compressed = KMerge::compress(sub);
    std::cout << "Finished compressing data" << std::endl;
    comp_size.push_back(compressed.size());

    std::cout << "Writing data for part " << i << std::endl;
    rc = unqlite_kv_store(db,(key1 + std::string("part|") + std::to_string(i)).c_str(),-1,&compressed[0],sizeof(uint)*comp_size[i]);
    REQUIRE(rc == UNQLITE_OK);
    std::cout << "Finished writing data" << std::endl;
    sub.clear();
  }
  std::vector<uint> mydataback;
  for (uint i=0; i < ceil((double) N/part_size);i++) {
    std::cout << "Reading data for " << i << std::endl;
    n_bytes = comp_size[i]*sizeof(uint);
    value = new uint [comp_size[i]];
    rc = unqlite_kv_fetch(db,(key1 + std::string("part|") + std::to_string(i)).c_str(),-1, value, &n_bytes);
    REQUIRE(rc == UNQLITE_OK);
    std::cout << "Finished reading data" << std::endl;
   
    std::cout << "Copying data" << std::endl;
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > compressed((uint*) &value[0], (uint*) &value[0] + n_bytes/sizeof(uint));
    std::cout << "Finished copying data" << std::endl;

    delete [] value;
    std::cout << "Uncompressing data" << std::endl;
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > sub = KMerge::uncompress(compressed, N);
    std::cout << "Finished uncompressing data" << std::endl;
    std::cout << "Appending data for " << i << std::endl;
    mydataback.insert(mydataback.end(), sub.begin(), sub.end());
    std::cout << "Finished appending data" << std::endl;
  }


  std::cout << "Checking data" << std::endl;
  REQUIRE(mydata.size() == mydataback.size());
  std::cout << "Finished checking data" << std::endl;

  unqlite_close(db);

  if( remove( "test.db" ) != 0 )
    perror( "Error deleting file" );

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

  std::cout << "Current memory usage in GB: " << KMerge::memory_used() << std::endl;
}


TEST_CASE("HashMapDump", "ConcurrentTest") {
  ulib::chain_hash_map<uint, uint> d_map(100000000);
  btree::btree_map<uint, uint> l_map;
  std::string filename("hashes.bin");
  
  d_map[0] = 0;
  d_map[1] = 1;
  d_map[2] = 2;
  d_map[3] = 3;

  KMerge::dump_hashes(d_map, filename);

  KMerge::load_hashes(l_map, filename);

  uint i = 0;
  for (btree::btree_map<uint, uint>::const_iterator iter = l_map.begin(); iter != l_map.end(); iter++) {
    REQUIRE(l_map[i] == i);
    i++;
  }

  if( remove( filename.c_str() ) != 0 )
    perror( "Error deleting file" );

}

TEST_CASE("QuickTest", "[Test]") {
  std::set<uint> t;

  for(uint i = 0; i < MAX_UINT_VAL; i++) {
    t.insert(i);
  }

  REQUIRE(t.size() == MAX_UINT_VAL);

  std::cout << "Current memory usage in GB: " << KMerge::memory_used() << std::endl;
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
  std::vector<uint> hashes_out, counts_out;
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >comp_counts;
  ulib::chain_hash_map<uint, uint> hashed_counts(1000000);
  const string kmer1("AAAAA");
  const string kmer2("TTTTT");
  const uint kmer1_count = 84;
  const uint kmer2_count = 13;
  uint kmer1_pos = 0;
  uint kmer2_pos = 0;
  uint pos = 0;
  param_struct params;


  params.k_val_start = 5;
  params.k_val_end = 5;
  params.db_filename = "/home/darryl/Development/kmerge/tests/parse_example.db";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/genome.test.fasta.gz";
  params.group_name = "org1";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
  params.priority = 1;
  params.lock_filename = params.db_filename + std::string(".lck");
  


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


  int rc = unqlite_open(&(params.db), params.db_filename.c_str(), UNQLITE_OPEN_CREATE);
  REQUIRE(rc == UNQLITE_OK);


  kmerge->add_dataset(hashes_out, params.group_name + std::string("|kmer_hash"), params.db, false);
  kmerge->add_dataset(counts_out, params.group_name + std::string("|count"), params.db, true);

  delete kmerge;

  unqlite_int64 n_bytes;
  n_bytes = sizeof(uint);
  uint size;
  rc = unqlite_kv_fetch(params.db,(params.group_name + std::string("|kmer_hash|size")).c_str(),-1, &size, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);

  uint *value = new uint[size];
  n_bytes = size*sizeof(uint);
  rc = unqlite_kv_fetch(params.db,(params.group_name + std::string("|kmer_hash")).c_str(),-1,value, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);
  std::vector<uint> hashes_in;
  hashes_in.assign((uint*) &value[0], (uint*) &value[0] + size);
  delete [] value;
 
  n_bytes = sizeof(uint);
  uint compressed_size;
  rc = unqlite_kv_fetch(params.db,(params.group_name + std::string("|count|size")).c_str(),-1, &compressed_size, &n_bytes);
  REQUIRE(rc == UNQLITE_OK);
  value = new uint[compressed_size];
  n_bytes = sizeof(uint)*compressed_size;
  rc = unqlite_kv_fetch(params.db,(params.group_name + std::string("|count")).c_str(),-1,value, &n_bytes);
  comp_counts.assign((uint*) &value[0], (uint*) &value[0] + compressed_size);
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > counts_in = KMerge::uncompress(comp_counts, size);
  comp_counts.clear();
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE>>().swap(comp_counts);
  delete [] value;

  unqlite_close(params.db);


  for (pos = 0; pos < hashes_in.size(); pos++) {
    REQUIRE(hashes_out[pos] == hashes_in[pos]);
    REQUIRE(counts_out[pos] == counts_in[pos]);
  }

  if( remove( params.db_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

}

TEST_CASE("ParseKmerCountsDumpMapAndCreateDB", "[HashTest]") {
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >comp_counts;
  ulib::chain_hash_map<uint, uint> hashed_counts(100000000);
  const string kmer1("AAAAA");
  const string kmer2("TTTTT");
  const uint kmer1_count = 84;
  const uint kmer2_count = 13;
  uint kmer1_pos = 0;
  uint kmer2_pos = 0;
  uint pos = 0;
  param_struct params;
  dlib::thread_pool tp(7);
  std::mutex mtx;
  std::condition_variable cv;

  params.k_val_start = 5;
  params.k_val_end = 13;
  params.db_filename = "/home/darryl/Development/kmerge/tests/dump_example.db";
  params.seq_filename = "/home/darryl/Development/kmerge/tests/208831/208831.fasta.gz";
  params.group_name = "208831";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
  params.priority = 1;
  params.lock_filename = params.db_filename + std::string(".lck");

  params.dump_mtx = &mtx;
  params.cv = &cv;
  params.ready = true;
  params.finished_hashing = false;
  params.hashed_counts = &hashed_counts;
  params.is_ref = true;
  params.dump_filename = "hashes.bin";
  params.writing.resize(params.num_threads, 0);
  params.polling_done = false;

  KMerge* kmerge = new KMerge(params.db_filename, "lookup3", ".", 0.5);
  params.kmerge = kmerge;
  
  KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
  tp.add_task(*task, &KMerge::BuilderTask::check_memory); 
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
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > counts_in = KMerge::uncompress(comp_counts, size);
  comp_counts.clear();
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE>>().swap(comp_counts);
  delete [] value;

  unqlite_close(params.db);

  REQUIRE(hashes_in.size() == 6091121);
  REQUIRE(counts_in.size() == 6091121);


  if( remove( params.db_filename.c_str() ) != 0 )
    perror( "Error deleting file" );
  if( remove( params.dump_filename.c_str() ) != 0 )
    perror( "Error deleting file" );
 

}


TEST_CASE("ThreadedParseKmerCountsAndCreateDBFromFastq", "[HashTest]") {
  param_struct params;
  std::string value;
  std::stringstream ss_delete;
  std::vector<uint> comp_hashes;
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > comp_counts;
  dlib::thread_pool tp(5);
  std::mutex mtx;
  std::condition_variable cv;
  ulib::chain_hash_map<uint, uint> hashed_counts(100000000);

  params.k_val_start = 3;
  params.k_val_end = 7;
  params.db_filename = "thread_fq.db";
  params.lock_filename = params.db_filename + std::string(".lck");
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz";
  params.group_name = "sample";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
  params.is_ref = false;
  params.dump_mtx = &mtx;
  params.cv = &cv;
  params.ready = true;
  params.finished_hashing = false;
  params.hashed_counts = &hashed_counts;
  params.dump_filename = "hashes.bin";
  params.writing.resize(params.num_threads, 0);
  params.polling_done = false;

  KMerge *kmerge = new KMerge(params.db_filename, "lookup3", ".");

  params.kmerge = kmerge;

  KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
  tp.add_task(*task, &KMerge::BuilderTask::check_memory);
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
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > counts_in;

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
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > part_counts = KMerge::uncompress(comp_counts, size);
    comp_counts.clear();
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >().swap(comp_counts);
    delete [] value;

    
    hashes_in.insert(hashes_in.end(), part_hashes.begin(), part_hashes.end());
    counts_in.insert(counts_in.end(), part_counts.begin(), part_counts.end());
  }
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
  if( remove( params.dump_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

}



TEST_CASE("ThreadedParseKmerCountsAndCreateDB", "[HashTest]") {
  std::vector<uint> hashes_in, comp_hashes;
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > counts_in, comp_counts;
  uint kmer1_count, kmer2_count, kmer1_pos, kmer2_pos, pos, total_uncompressed_size;
  param_struct params1, params2, params3;
  std::map<std::string, std::string> taxonomy;
  std::string line;
  char* value;
  ifstream in_file;
  std::stringstream ss_in, ss_delete;
  unqlite_int64 n_bytes;
  uint num_parts;
  ulib::chain_hash_map<uint, uint> hashed_counts1(100000000), hashed_counts2(100000000), hashed_counts3(100000000);

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
  params1.dump_filename = "hashes1.bin";
  params1.num_threads = (params1.k_val_end - params1.k_val_start) / 2 + 1;
  params1.priority = 1;
  params1.is_ref = true;
  params1.dump_mtx = new std::mutex;
  params1.cv = new std::condition_variable;
  params1.ready = true;
  params1.finished_hashing = false;
  params1.hashed_counts = &hashed_counts1;
  params1.writing.resize(params1.num_threads, 0);
  params1.polling_done = false;

  params2.seq_filename = "/home/darryl/Development/kmerge/tests/209328/209328.fasta.gz";
  params2.group_name = "209328";
  params2.dump_filename = "hashes2.bin";
  params2.num_threads = (params2.k_val_end - params2.k_val_start) / 2 + 1;
  params2.priority = 2;
  params2.is_ref = true;
  params2.dump_mtx = new std::mutex;
  params2.cv = new std::condition_variable;
  params2.ready = true;
  params2.finished_hashing = false;
  params2.hashed_counts = &hashed_counts2;
  params2.writing.resize(params2.num_threads, 0);
  params2.polling_done = false;

  params3.seq_filename = "/home/darryl/Development/kmerge/tests/54095/54095.fasta.gz";
  params3.group_name = "54095";
  params3.dump_filename = "hashes3.bin";
  params3.num_threads = (params3.k_val_end - params3.k_val_start) / 2 + 1;
  params3.priority = 3;
  params3.is_ref = true;
  params3.dump_mtx = new std::mutex;
  params3.cv = new std::condition_variable;
  params3.ready = true;
  params3.finished_hashing = false;
  params3.hashed_counts = &hashed_counts3;
  params3.writing.resize(params3.num_threads, 0);
  params3.polling_done = false;

  uint thread_count = 6;

  dlib::thread_pool tp(thread_count);



  KMerge* kmerge = new KMerge(params1.db_filename, "lookup3", ".", 0.5);

  params1.kmerge = kmerge;
  params2.kmerge = kmerge;
  params3.kmerge = kmerge;

  KMerge::BuilderTask t1(params1);
  KMerge::BuilderTask t2(params2);
  KMerge::BuilderTask t3(params3);

  tp.add_task(t1, &KMerge::BuilderTask::check_memory);
  tp.add_task(t1, &KMerge::BuilderTask::execute);
  tp.add_task(t2, &KMerge::BuilderTask::check_memory);
  tp.add_task(t2, &KMerge::BuilderTask::execute);
  tp.add_task(t3, &KMerge::BuilderTask::check_memory);
  tp.add_task(t3, &KMerge::BuilderTask::execute);

  tp.wait_for_all_tasks();

  delete kmerge;
  delete params1.dump_mtx;
  delete params1.cv;
  delete params2.dump_mtx;
  delete params2.cv;
  delete params3.dump_mtx;
  delete params3.cv;

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
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > part_counts = KMerge::uncompress(comp_counts, size);
    comp_counts.clear();
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >().swap(comp_counts);
    delete [] value;

    
    hashes_in.insert(hashes_in.end(), part_hashes.begin(), part_hashes.end());
    counts_in.insert(counts_in.end(), part_counts.begin(), part_counts.end());
    
  }

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
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >().swap(counts_in);
 
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
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > part_counts = KMerge::uncompress(comp_counts, size);
    comp_counts.clear();
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >().swap(comp_counts);
    delete [] value;

    
    hashes_in.insert(hashes_in.end(), part_hashes.begin(), part_hashes.end());
    counts_in.insert(counts_in.end(), part_counts.begin(), part_counts.end());
    
  }



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
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >().swap(counts_in);

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
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > part_counts = KMerge::uncompress(comp_counts, size);
    comp_counts.clear();
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >().swap(comp_counts);
    delete [] value;

    
    hashes_in.insert(hashes_in.end(), part_hashes.begin(), part_hashes.end());
    counts_in.insert(counts_in.end(), part_counts.begin(), part_counts.end());
  }


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
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >().swap(counts_in);

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
  if( remove( params1.dump_filename.c_str() ) != 0 )
    perror( "Error deleting file" );
  if( remove( params2.dump_filename.c_str() ) != 0 )
    perror( "Error deleting file" );
  if( remove( params3.dump_filename.c_str() ) != 0 )
    perror( "Error deleting file" );


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
  std::string line;
  char* value;
  std::stringstream ss_in;
  std::map<std::string, std::string> taxonomy;

  unqlite *db;

  int rc = unqlite_open(&db, db_filename.c_str(), UNQLITE_OPEN_CREATE);
  REQUIRE(rc == UNQLITE_OK);

  
  KMerge *kmerge = new KMerge(db_filename, "lookup3", ".");

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

