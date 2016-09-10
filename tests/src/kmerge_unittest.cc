#define CATCH_CONFIG_MAIN

#include <limits.h>
#include "kmerge.h"
#include <math.h>
#include "catch.hpp"
#include <sstream>
#include <dlib/threads.h>
#include <dlib/serialize.h>
#include <dlib/logger.h>
#include <sys/stat.h>
#include <iomanip>


TEST_CASE("BasicStringKernelTest", "[SSKernel]") {
  unsigned long long x, y;
  int i, k=5;
  std::vector<std::string> data;
  std::string kmer;
  uint count = 0;

  // only generate the lexicographically ordered k-mers to reduce the total number 
  // needed in the kernel by 2

  for (x = 0; x < 1ULL<<(2*k); ++x) {
    for (i = 0, y = x; i < k; ++i, y >>= 2)
      kmer.push_back("acgt"[y&3]);
    count++;
    // get reverse complement of k-mer
    std::string rc = KMerge::rev_comp(kmer);
    // check which is lexicographically first
    // if this is lexicographically first, add k-mer, if not do nothing
    if (kmer < rc) data.push_back(kmer);
    kmer.clear();
  }
 
  // Kernel parameters 
  const float c = 1e12;
  const int normalize = 1;
  const int symbol_size = 255;  // A size of an alphabet
  const int max_length = 1000;  // A maximum sequence length
  int kn = 3;                   // A level of susbsequence matching
  double lambda = 0.5;          // A decay factor

  string_kernel sk(c, normalize, symbol_size, max_length, kn, lambda);
  sk.set_data(data);
  sk.compute_kernel();

  REQUIRE(sk(std::string("AAAAA"), std::string("aaaaa")) > 0);
  REQUIRE(sk(std::string("AAAAA"), std::string("AAAAA")) > sk(std::string("AAAAA"), std::string("AACCA")));
  REQUIRE(sk(std::string("AAAAA"), std::string("AAAAA")) == sk(std::string("AAAAA"), std::string("TTTTT")));
    
}

TEST_CASE("SeralizeSSKernel", "[SSKernel]") {
  unsigned long long x, y;
  int i, k=5;
  std::vector<std::string> data;
  std::string kmer;
  uint count = 0;
  std::ofstream ofile("kernel.bin");

  // only generate the lexicographically ordered k-mers to reduce the total number 
  // needed in the kernel by 2

  for (x = 0; x < 1ULL<<(2*k); ++x) {
    for (i = 0, y = x; i < k; ++i, y >>= 2)
      kmer.push_back("acgt"[y&3]);
    count++;
    // get reverse complement of k-mer
    std::string rc = KMerge::rev_comp(kmer);
    // check which is lexicographically first
    // if this is lexicographically first, add k-mer, if not do nothing
    if (kmer < rc) data.push_back(kmer);
    kmer.clear();
  }
 
  // Kernel parameters 
  const float c = 1e12;
  const int normalize = 1;
  const int symbol_size = 255;  // A size of an alphabet
  const int max_length = 1000;  // A maximum sequence length
  int kn = 2;                   // A level of susbsequence matching
  double lambda = 0.5;          // A decay factor

  string_kernel sk(c, normalize, symbol_size, max_length, kn, lambda);
  sk.set_data(data);
  sk.compute_kernel();
  serialize(sk, ofile);

  std::ifstream ifile("kernel.bin");
  string_kernel sk2;
  deserialize(sk2, ifile);
}



TEST_CASE("GenerateTruncatedHashes", "[KmerGenerator]") {
  std::set<uint> murmur_set;
  int i, k=7;
  unsigned long long x, y;
  uint bits = 8;
  
  std::string kmer;
  uint count = 0;
  
  for (x = 0; x < 1ULL<<(2*k); ++x) {
    for (i = 0, y = x; i < k; ++i, y >>= 2)
      kmer.push_back("ACGT"[y&3]);
    count++;
    murmur_set.insert((uint8_t) KMerge::hash_kmer(kmer, MURMUR, BITS_8));
    kmer.clear();
  }
  
  REQUIRE(pow(4,k) == count);
  REQUIRE(pow(2,bits) >= murmur_set.size());

  murmur_set.clear();
  k=9;
  bits = 16;
  count = 0;

  for (x = 0; x < 1ULL<<(2*k); ++x) {
    for (i = 0, y = x; i < k; ++i, y >>= 2)
      kmer.push_back("ACGT"[y&3]);
    count++;
    murmur_set.insert((uint16_t) KMerge::hash_kmer(kmer, MURMUR, BITS_16));
    kmer.clear();
  }

  REQUIRE(pow(4,k) == count);
  REQUIRE(pow(2,bits) >= murmur_set.size());

}



TEST_CASE("GenerateAllKMersOfGivenLength", "[KmerGenerator]") {
  std::set<uint> lookup3_set, city_set, murmur_set;
  int i, k=7;
  unsigned long long x, y;
  
  std::string kmer;
  uint count = 0;
  
  for (x = 0; x < 1ULL<<(2*k); ++x) {
    for (i = 0, y = x; i < k; ++i, y >>= 2)
      kmer.push_back("ACGT"[y&3]);
    count++;
    lookup3_set.insert(KMerge::hash_kmer(kmer, LOOKUP3, BITS_32));
    city_set.insert(KMerge::hash_kmer(kmer, MURMUR, BITS_32));
    murmur_set.insert(KMerge::hash_kmer(kmer, CITY, BITS_32));
    kmer.clear();
  }
  
  REQUIRE(pow(4,k) == count);
  REQUIRE((pow(4,k) / 2) == lookup3_set.size());
  REQUIRE((pow(4,k) / 2) == city_set.size());
  REQUIRE((pow(4,k) / 2) == murmur_set.size());

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

TEST_CASE("CountKmersInPEFastqFile", "[HashTest]") {

  int k = 5, l;
  std::vector<std::map<std::string,uint> > kmer_counts;
  std::map<std::string, uint> jf_counts;
  kseq_t *seq;
  gzFile fp;
  std::vector<std::string> kmers;

  fp = gzopen("/home/darryl/Development/kmerge/tests/sample/sample.pe.fastq.gz", "rb");
  seq = kseq_init(fp);
  uint seq_id = 0;
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
    if (seq_id % 2 == 1) {
      std::map<std::string, uint> cnt_map;
      for (auto it = kmers.begin(); it != kmers.end(); it++) {
	cnt_map[*it]++;
      }
      kmer_counts.push_back(cnt_map);
      kmers.clear();
      cnt_map.clear();
    }
    seq_id++;
  }
  kseq_destroy(seq);
  gzclose(fp);

  std::string base_dir("/home/darryl/Development/kmerge/tests/sample/split/");
  ifstream in_file;
  for (uint i=0; i<kmer_counts.size(); i++) {
    std::stringstream filename;
    filename << "x" << setfill('0') << setw(4) << i << ".cnt";
    in_file.open(base_dir + filename.str());
    std::string seq2;
    uint count;
    if (in_file.is_open()) {
      while ( !in_file.eof() ) {
	in_file >> seq2 >> count;
	if (seq2 != "") {
	  jf_counts[seq2] = count;
	  // multiply count by 2 because test PE reads are simply duplicates
	  REQUIRE(kmer_counts[i][seq2] == count*2);
	}
      }
    }    
    in_file.close();
 
    REQUIRE(jf_counts.size() == kmer_counts[i].size());
    jf_counts.clear();
  }
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
  params.seq_filename = "/home/darryl/Development/kmerge/tests/208831/208831.fasta.gz";
  params.group_name = "208831";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;  



  //~/bin/kanalyze-0.9.5/count -d 20 -f fastagz -k 3 -o 208831/sample.k3.txt 208831/208831.fasta.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastagz -k 5 -o 208831/sample.k5.txt 208831/208831.fasta.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastagz -k 7 -o 208831/sample.k7.txt 208831/208831.fasta.gz

  files.push_back("/home/darryl/Development/kmerge/tests/208831/sample.k3.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/208831/sample.k5.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/208831/sample.k7.txt");

  KMerge *kmerge = new KMerge("lookup3", ".", "./reference", 0);

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
  btree::btree_map<std::string, btree::btree_map<uint,uint> > hashed_counts;
  btree::btree_map<uint, uint> combined_counts, jf_hashed_counts;
  param_struct params;
  std::vector<std::string> files;
  kseq_t *seq;
  gzFile fp;
  int l;
  std::mutex mtx;


  params.k_val_start = 3;
  params.k_val_end = 7;
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz";
  params.group_name = "sample";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;  


  //~/bin/kanalyze-0.9.5/count -d 20 -f fastqgz -k 3 -o sample/sample.k3.txt sample/sample.fastq.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastqgz -k 5 -o sample/sample.k5.txt sample/sample.fastq.gz
  //~/bin/kanalyze-0.9.5/count -d 20 -f fastqgz -k 7 -o sample/sample.k7.txt sample/sample.fastq.gz

  files.push_back("/home/darryl/Development/kmerge/tests/sample/sample.k3.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/sample/sample.k5.txt");
  files.push_back("/home/darryl/Development/kmerge/tests/sample/sample.k7.txt");


  KMerge *kmerge = new KMerge("lookup3", ".", "./reference", 0);

  params.kmerge = kmerge;

  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  uint seq_idx = 0;
  std::vector<std::string> seq_tup;
  std::string base_id;
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s), seq_name(seq->name.s);
    seq_tup.push_back(seq_str);
    if (seq_idx % 1 == 0) {
      for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
	bool success = params.kmerge->hash_seq(seq_tup, k, hashed_counts, seq_name, mtx);
	REQUIRE(success == true);
      }
      seq_tup.clear();
    }
    seq_idx++;
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

  for (auto s_iter: hashed_counts) {
    for (auto m_iter: s_iter.second) {
      combined_counts[m_iter.first] += m_iter.second;
    }
  }

  REQUIRE(jf_hashed_counts.size() == combined_counts.size());

  for (btree::btree_map<uint, uint>::const_iterator b_it = jf_hashed_counts.begin(); b_it != jf_hashed_counts.end(); b_it++) {
    REQUIRE(combined_counts[b_it->first] == b_it->second);
  }

}


TEST_CASE("GeneratePairedEndBaseID", "[HashTest]") {
  std::string r_seq = R"(sequence\1)", l_seq = R"(sequence\2)";
  std::string base;

  base = KMerge::get_seq_base_id(r_seq, l_seq);
  REQUIRE(base == "sequence");

  r_seq = R"(sequence|1)";
  l_seq = R"(sequence|2)";
  base = KMerge::get_seq_base_id(r_seq, l_seq);
  REQUIRE(base == "sequence");

  r_seq = R"(sequence:1)";
  l_seq = R"(sequence:2)";
  base = KMerge::get_seq_base_id(r_seq, l_seq);
  REQUIRE(base == "sequence");

}

TEST_CASE("CountHashedKmersInPEFastqFile", "[HashTest]") {
  btree::btree_map<std::string, btree::btree_map<uint,uint> > hashed_counts;
  btree::btree_map<uint, uint> jf_hashed_counts;
  param_struct params;
  std::vector<std::string> seq_ids;
  kseq_t *seq;
  gzFile fp;
  int l;
  std::mutex mtx;


  params.k_val_start = 5;
  params.k_val_end = 5;
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.pe.fastq.gz";
  params.group_name = "sample";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;  



  KMerge *kmerge = new KMerge("lookup3", ".", "./reference", 0);

  params.kmerge = kmerge;
 
  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  uint seq_idx = 0;
  std::vector<std::string> seq_tup, seq_names;
  std::string base_id;
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s), seq_name(seq->name.s);
    seq_tup.push_back(seq_str);
    seq_names.push_back(seq_name);
    if (seq_idx % 2 == 1) {
      base_id = KMerge::get_seq_base_id(seq_names[0], seq_names[1]);
      REQUIRE(base_id == std::string("sequence") + std::to_string(hashed_counts.size()+1));
      seq_ids.push_back(base_id);
      for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
	bool success = params.kmerge->hash_seq(seq_tup, k, hashed_counts, base_id, mtx);
	REQUIRE(success == true);
      }
      seq_tup.clear();
      seq_names.clear();
    }
    seq_idx++;
  }
  kseq_destroy(seq);
  gzclose(fp);


  std::string base_dir("/home/darryl/Development/kmerge/tests/sample/split/");
  std::string seq_str2, last = "";
  uint count;
  ifstream in_file;
  
  for (auto seq_id: seq_ids) {
    std::vector<std::string>::iterator iter = std::find(seq_ids.begin(), seq_ids.end(), seq_id);
    uint idx = iter - seq_ids.begin();
    std::stringstream filename;
    filename << "x" << setfill('0') << setw(4) << idx << ".cnt";
    in_file.open(base_dir + filename.str());
    std::string seq2;
    uint count;
 
    if (in_file.is_open()) {
      while ( !in_file.eof() ) {
	in_file >> seq_str2 >> count;
	if (seq_str2 != last) {
	  uint hash = kmerge->hash_kmer(seq_str2);
	  jf_hashed_counts[hash] = jf_hashed_counts[hash] + count;
	}
	last = seq_str2;
      }
    }
    in_file.close();
    REQUIRE(jf_hashed_counts.size() == hashed_counts[seq_id].size());
    for (btree::btree_map<uint, uint>::const_iterator b_it = jf_hashed_counts.begin(); b_it != jf_hashed_counts.end(); b_it++) {
      REQUIRE(hashed_counts[seq_id][b_it->first] == b_it->second*2);
    }
    jf_hashed_counts.clear();
  }

  delete kmerge;

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
  params.seq_filename = "/zenodotus/masonlab/pathomap_scratch/darryl/k-mer/data/208831/208831.fasta.gz";
  params.group_name = "208831";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  KMerge *kmerge = new KMerge("lookup3", ".", "./reference", 0);
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
  btree::btree_map<std::string, btree::btree_map<uint, uint> > hashed_counts, hashed_counts2;
  int l;
  param_struct params;
  std::vector<std::tuple<std::vector<std::string>, std::string, uint> > jobs;
  std::mutex mtx;
  kseq_t *seq;
  gzFile fp;
  std::vector<std::string> seq_ids, seq_tup;
  uint seq_idx = 0;
  
  params.k_val_start = 3;
  params.k_val_end = 11;
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.fastq.gz";
  params.group_name = "sample";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  KMerge *kmerge = new KMerge("lookup3", ".", "./reference", 0);
  params.kmerge = kmerge;


  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s), seq_name(seq->name.s);
    seq_tup.push_back(seq_str);
    seq_ids.push_back(seq_name);
    if (seq_idx % 1 == 0) {
      for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
	jobs.push_back(std::make_tuple(seq_tup, seq_name, k));
      }
      seq_tup.clear();
    }
    seq_idx++;
  }
  REQUIRE(jobs.size() == 5000);
  KMerge::HashSeqs func(params, hashed_counts, jobs, mtx, false);
  dlib::parallel_for(params.num_threads, 0, jobs.size(), func);
  jobs.clear();
  kseq_destroy(seq);
  gzclose(fp);
  
 
  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  seq_idx = 0;
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s), seq_name(seq->name.s);
    seq_tup.push_back(seq_str);
    if (seq_idx % 1 == 0) {
      for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
	bool success = params.kmerge->hash_seq(seq_tup, k, hashed_counts2, seq_name, mtx);
	REQUIRE(success == true);
      }
      seq_tup.clear();
    }
    seq_idx++;
  }
  kseq_destroy(seq);
  gzclose(fp);


  REQUIRE(hashed_counts.size() == hashed_counts2.size());

  for (auto seq_id: seq_ids) {
    REQUIRE(hashed_counts[seq_id].size() == hashed_counts2[seq_id].size());
    for (auto cnt_map: hashed_counts[seq_id]) {
      uint hash = cnt_map.first;
      uint count = cnt_map.second;
      REQUIRE(hashed_counts2[seq_id][hash] == count);
    }

  }

  delete kmerge;
}


TEST_CASE("CountHashedKmersInParallelPEFastq", "[HashTest]") {
  btree::btree_map<std::string, btree::btree_map<uint, uint> > hashed_counts, hashed_counts2;
  int l;
  param_struct params;
  std::vector<std::tuple<std::vector<std::string>, std::string, uint> > jobs;
  std::mutex mtx;
  kseq_t *seq;
  gzFile fp;
  std::vector<std::string> seq_tup, seq_ids, seq_names;
  uint seq_idx = 0;
  std::string base_id;

  
  params.k_val_start = 5;
  params.k_val_end = 5;
  params.seq_filename = "/home/darryl/Development/kmerge/tests/sample/sample.pe.fastq.gz";
  params.group_name = "sample";
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;

  KMerge *kmerge = new KMerge("lookup3", ".", "./reference", 0);
  params.kmerge = kmerge;


  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s), seq_name(seq->name.s);
    seq_tup.push_back(seq_str);
    seq_names.push_back(seq_name);
    if (seq_idx % 2 == 1) {
      base_id = KMerge::get_seq_base_id(seq_names[0], seq_names[1]);
      seq_ids.push_back(base_id);
      REQUIRE(base_id == std::string("sequence") + std::to_string(seq_ids.size()));
      for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
	jobs.push_back(std::make_tuple(seq_tup, base_id, k));
      }
      seq_tup.clear();
      seq_names.clear();
    }
    seq_idx++;
  }
  KMerge::HashSeqs func(params, hashed_counts, jobs, mtx, false);
  dlib::parallel_for(params.num_threads, 0, jobs.size(), func);
  jobs.clear();
  kseq_destroy(seq);
  gzclose(fp);
  

  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  seq_idx = 0;
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s), seq_name(seq->name.s);
    seq_tup.push_back(seq_str);
    seq_names.push_back(seq_name);
    if (seq_idx % 2 == 1) {
      base_id = KMerge::get_seq_base_id(seq_names[0], seq_names[1]);
      REQUIRE(base_id == std::string("sequence") + std::to_string(hashed_counts2.size()+1));
      for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
	bool success = params.kmerge->hash_seq(seq_tup, k, hashed_counts2, base_id, mtx);
	REQUIRE(success == true);
      }
      seq_tup.clear();
      seq_names.clear();
    }
    seq_idx++;
  }
  kseq_destroy(seq);
  gzclose(fp);


  REQUIRE(hashed_counts.size() == hashed_counts2.size());

  for (auto seq_id: seq_ids) {
    REQUIRE(hashed_counts[seq_id].size() == hashed_counts2[seq_id].size());
    for (auto cnt_map: hashed_counts[seq_id]) {
      uint hash = cnt_map.first;
      uint count = cnt_map.second;
      REQUIRE(hashed_counts2[seq_id][hash] == count);
    }

  }

  delete kmerge;
}



TEST_CASE("ParseKmerCountsAndCreateDB", "[HashTest]") {
  std::vector<uint> hashes_in, counts_in;
  const string kmer1("AAAAA");
  const string kmer2("TTTTT");
  const uint kmer1_count = 84;
  const uint kmer2_count = 13;
  uint kmer1_pos = 0;
  uint kmer2_pos = 0;
  uint pos = 0;
  param_struct params;
  dlib::thread_pool tp(1);
  cs decompressor;
  std::stringstream ss;
  std::ifstream fs;

  params.k_val_start = 5;
  params.k_val_end = 13;
  params.group_name = "208831";
  params.seq_filename = std::string("./") + params.group_name + std::string("/") + params.group_name + std::string(".fasta.gz");
  params.hashes_filename = std::string("./reference/") + params.group_name + std::string(".hashes.bin");
  params.counts_filename = std::string("./reference/") + params.group_name + std::string(".counts.bin");
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
  params.is_ref = true;

  KMerge* kmerge = new KMerge("lookup3", ".", "./reference", 0);
  params.kmerge = kmerge;
  
  KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
  tp.add_task(*task, &KMerge::BuilderTask::execute);

  tp.wait_for_all_tasks();
 
  delete kmerge;

  fs.open (params.hashes_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(hashes_in, ss);
  ss.str(std::string());
  std::partial_sum(hashes_in.begin(), hashes_in.end(), hashes_in.begin());

  fs.open (params.counts_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(counts_in, ss);
  ss.str(std::string());

  REQUIRE(hashes_in.size() == 6091121);
  REQUIRE(counts_in.size() == 6091121);


  if( remove( params.hashes_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

  if( remove( params.counts_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

  if( remove( std::string("./reference/" + params.group_name + ".taxonomy.bin").c_str() ) != 0 )
    perror( "Error deleting file" );
}


TEST_CASE("ThreadedParseKmerCountsAndCreateDBFromFastq", "[HashTest]") {
  param_struct params;
  std::vector<uint> hashes_in, counts_in, indices_in;
  std::vector<std::string> ids_in, seq_tup;
  btree::btree_map<std::string, btree::btree_map<uint,uint> > hashed_counts;
  dlib::thread_pool tp(1);
  std::mutex mtx;
  cs decompressor;
  std::stringstream ss;
  std::ifstream fs;
  int l;
  kseq_t *seq;
  gzFile fp;


  params.k_val_start = 3;
  params.k_val_end = 7;
  params.group_name = "sample";
  params.seq_filename = std::string("./") + params.group_name + std::string("/") + params.group_name + std::string(".fastq.gz");
  params.hashes_filename = std::string("./reference/") + params.group_name + std::string(".hashes.bin");
  params.counts_filename = std::string("./reference/") + params.group_name + std::string(".counts.bin");
  params.indices_filename = std::string("./reference/") + params.group_name + std::string(".indices.bin");
  params.ids_filename = std::string("./reference/") + params.group_name + std::string(".ids.bin");
  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
  params.is_ref = false;
  params.paired_end = false;

  KMerge *kmerge = new KMerge("lookup3", ".", "./reference", 0);

  params.kmerge = kmerge;

  KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
  tp.add_task(*task, &KMerge::BuilderTask::execute);
  tp.wait_for_all_tasks();

  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  uint seq_idx = 0;
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s), seq_name(seq->name.s);
    seq_tup.push_back(seq_str);
    if (seq_idx % 1 == 0) {
      for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
	bool success = params.kmerge->hash_seq(seq_tup, k, hashed_counts, seq_name, mtx);
	REQUIRE(success == true);
      }
      seq_tup.clear();
    }
    seq_idx++;
  }
  kseq_destroy(seq);
  gzclose(fp);


  delete kmerge;

  fs.open (params.hashes_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(hashes_in, ss);
  ss.str(std::string());

  fs.open (params.counts_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(counts_in, ss);
  ss.str(std::string());

  fs.open (params.indices_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(indices_in, ss);
  ss.str(std::string());

  for (uint i=0; i < indices_in.size() - 1; i++) {
    std::partial_sum(hashes_in.begin() + indices_in[i], hashes_in.begin() + indices_in[i+1], hashes_in.begin() + indices_in[i]);
  }

  fs.open (params.ids_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(ids_in, ss);
  ss.str(std::string());


  REQUIRE(counts_in.size() == 289211);
  REQUIRE(hashes_in.size() == 289211);
  REQUIRE(indices_in.size() == 1001);
  REQUIRE(ids_in.size() == 1000);
  REQUIRE(indices_in.back() == 289211);
  uint total_hash_count=0;
  for (auto iter: hashed_counts) {
    total_hash_count += iter.second.size();
  }
  REQUIRE(total_hash_count == 289211);
 

  //ensure that hashes are in sorted order in each sample
  for (uint i=0; i < ids_in.size(); i++) {
    std::string seq_id = ids_in[i];
    uint last = 0;
    REQUIRE((indices_in[i+1]-indices_in[i]) == hashed_counts[seq_id].size());
    for (uint j=indices_in[i]; j < indices_in[i+1]; j++) {
      if (j != indices_in[i]) {
	REQUIRE(hashes_in[j] > last);
      }
      REQUIRE(hashed_counts[seq_id][hashes_in[j]] == counts_in[j]);
      last = hashes_in[j];
    }
  }
  

  if( remove( params.hashes_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

  if( remove( params.counts_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

  if( remove( params.indices_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

  if( remove( params.ids_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

}



TEST_CASE("ThreadedParseKmerCountsAndCreateDB", "[HashTest]") {
  std::vector<uint> hashes_in, counts_in;
  uint kmer1_count, kmer2_count, kmer1_pos, kmer2_pos, pos;
  param_struct params1, params2, params3;
  std::map<std::string, std::string> taxonomy;
  std::string line;
  ifstream in_file;
  cs decompressor;
  std::stringstream ss;
  std::ifstream fs;


  const string kmer1("AAAAA");
  const string kmer2("GCGAT");


  params1.k_val_start = 5;
  params1.k_val_end = 5;
  params2.k_val_start = 5;
  params2.k_val_end = 5;
  params3.k_val_start = 5;
  params3.k_val_end = 5;


  params1.group_name = "208831";
  params1.seq_filename = std::string("./") + params1.group_name + std::string("/") + params1.group_name + std::string(".fasta.gz");
  params1.hashes_filename = std::string("./reference/") + params1.group_name + std::string(".hashes.bin");
  params1.counts_filename = std::string("./reference/") + params1.group_name + std::string(".counts.bin");
  params1.num_threads = (params1.k_val_end - params1.k_val_start) / 2 + 1;
  params1.is_ref = true;

  params2.group_name = "209328";
  params2.seq_filename = std::string("./") + params2.group_name + std::string("/") + params2.group_name + std::string(".fasta.gz");
  params2.hashes_filename = std::string("./reference/") + params2.group_name + std::string(".hashes.bin");
  params2.counts_filename = std::string("./reference/") + params2.group_name + std::string(".counts.bin");
  params2.num_threads = (params2.k_val_end - params2.k_val_start) / 2 + 1;
  params2.is_ref = true;

  params3.group_name = "54095";
  params3.seq_filename = std::string("./") + params3.group_name + std::string("/") + params3.group_name + std::string(".fasta.gz");
  params3.hashes_filename = std::string("./reference/") + params3.group_name + std::string(".hashes.bin");
  params3.counts_filename = std::string("./reference/") + params3.group_name + std::string(".counts.bin");
 
  params3.num_threads = (params3.k_val_end - params3.k_val_start) / 2 + 1;
  params3.is_ref = true;

  uint thread_count = 3;

  dlib::thread_pool tp(thread_count);



  KMerge* kmerge = new KMerge("lookup3", ".", "./reference", 0);

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

  fs.open (params1.hashes_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(hashes_in, ss);
  ss.str(std::string());
  std::partial_sum(hashes_in.begin(), hashes_in.end(), hashes_in.begin());

  fs.open (params1.counts_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(counts_in, ss);
  ss.str(std::string());


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

  kmer1_count = 6150 /*AAAAA*/ + 6021 /*TTTTT*/;
  kmer2_count = 10775 /*GCGAT*/ + 10855 /*ATCGC*/;

  for (pos = 0; pos < 512 /*3^5*/; pos++) {
    if (hashes_in[pos] == KMerge::hash_kmer(kmer1, LOOKUP3, BITS_32)) {
      kmer1_pos = pos;
    }
    if (hashes_in[pos] == KMerge::hash_kmer(kmer2, LOOKUP3, BITS_32)) {
      kmer2_pos = pos;
    }
  }


  REQUIRE(hashes_in[kmer1_pos] == KMerge::hash_kmer(kmer1, LOOKUP3, BITS_32));
  REQUIRE(hashes_in[kmer2_pos] == KMerge::hash_kmer(kmer2, LOOKUP3, BITS_32));
  REQUIRE(counts_in[kmer1_pos] == kmer1_count);
  REQUIRE(counts_in[kmer2_pos] == kmer2_count);

  hashes_in.clear();
  std::vector<uint>().swap(hashes_in);
  counts_in.clear();
  std::vector<uint>().swap(counts_in);
 
  fs.open (std::string("./reference/" + params1.group_name + ".taxonomy.bin").c_str(), ios::binary);
  dlib::deserialize(taxonomy, fs);
  fs.close();

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


  fs.open (params2.hashes_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(hashes_in, ss);
  ss.str(std::string());
  std::partial_sum(hashes_in.begin(), hashes_in.end(), hashes_in.begin());

  fs.open (params2.counts_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(counts_in, ss);
  ss.str(std::string());


  //ensure that hashes are in sorted order
  for (std::vector<uint>::const_iterator v_iter = hashes_in.begin(); v_iter != hashes_in.end(); v_iter++) {
    if (v_iter != hashes_in.begin()) {
      REQUIRE(*v_iter > last);
    }
    last = *v_iter;
  }

  for (pos = 0; pos < 512 /*3^5*/; pos++) {
    if (hashes_in[pos] == KMerge::hash_kmer(kmer1, LOOKUP3, BITS_32)) {
      kmer1_pos = pos;
    }
    if (hashes_in[pos] == KMerge::hash_kmer(kmer2, LOOKUP3, BITS_32)) {
      kmer2_pos = pos;
    }
  }

  REQUIRE(hashes_in[kmer1_pos] == KMerge::hash_kmer(kmer1, LOOKUP3, BITS_32));
  REQUIRE(hashes_in[kmer2_pos] == KMerge::hash_kmer(kmer2, LOOKUP3, BITS_32));
  REQUIRE(counts_in[kmer1_pos] == kmer1_count);
  REQUIRE(counts_in[kmer2_pos] == kmer2_count);

  hashes_in.clear();
  std::vector<uint>().swap(hashes_in);
  counts_in.clear();
  std::vector<uint>().swap(counts_in);

  fs.open (std::string("./reference/" + params2.group_name + ".taxonomy.bin").c_str(), ios::binary);
  dlib::deserialize(taxonomy, fs);
  fs.close();
 
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

  fs.open (params3.hashes_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(hashes_in, ss);
  ss.str(std::string());
  std::partial_sum(hashes_in.begin(), hashes_in.end(), hashes_in.begin());

  fs.open (params3.counts_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(counts_in, ss);
  ss.str(std::string());


  //ensure that hashes are in sorted order
  for (std::vector<uint>::const_iterator v_iter = hashes_in.begin(); v_iter != hashes_in.end(); v_iter++) {
    if (v_iter != hashes_in.begin()) {
      REQUIRE(*v_iter > last);
    }
    last = *v_iter;
  }


  for (pos = 0; pos < 512 /*3^5*/; pos++) {
    if (hashes_in[pos] == KMerge::hash_kmer(kmer1, LOOKUP3, BITS_32)) {
      kmer1_pos = pos;
    }
    if (hashes_in[pos] == KMerge::hash_kmer(kmer2, LOOKUP3, BITS_32)) {
      kmer2_pos = pos;
    }
  }

  REQUIRE(hashes_in[kmer1_pos] == KMerge::hash_kmer(kmer1, LOOKUP3, BITS_32));
  REQUIRE(hashes_in[kmer2_pos] == KMerge::hash_kmer(kmer2, LOOKUP3, BITS_32));
  REQUIRE(counts_in[kmer1_pos] == kmer1_count);
  REQUIRE(counts_in[kmer2_pos] == kmer2_count);

  hashes_in.clear();
  std::vector<uint>().swap(hashes_in);
  counts_in.clear();
  std::vector<uint>().swap(counts_in);


  fs.open (std::string("./reference/" + params3.group_name + ".taxonomy.bin").c_str(), ios::binary);
  dlib::deserialize(taxonomy, fs);
  fs.close();


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

  if( remove( params1.hashes_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

  if( remove( params1.counts_filename.c_str() ) != 0 )
    perror( "Error deleting file" );
  
  if( remove( std::string("./reference/" + params1.group_name + ".taxonomy.bin").c_str() ) != 0)
    perror ("Error deleting file" );

  if( remove( params2.hashes_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

  if( remove( params2.counts_filename.c_str() ) != 0 )
    perror( "Error deleting file" );
  
  if( remove( std::string("./reference/" + params2.group_name + ".taxonomy.bin").c_str() ) != 0)
    perror ("Error deleting file" );

  if( remove( params3.hashes_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

  if( remove( params3.counts_filename.c_str() ) != 0 )
    perror( "Error deleting file" );
  
  if( remove( std::string("./reference/" + params3.group_name + ".taxonomy.bin").c_str() ) != 0)
    perror ("Error deleting file" );


}

TEST_CASE("ThreadedParseKmerCountsAndCreateDBFromPEFastq", "[HashTest]") {
  param_struct params;
  btree::btree_map<std::string, btree::btree_map<uint,uint> > hashed_counts;
  std::vector<uint> hashes_in, counts_in, indices_in;
  std::vector<std::string> ids_in;
  dlib::thread_pool tp(1);
  std::mutex mtx;
  cs decompressor;
  std::stringstream ss;
  std::ifstream fs;
  std::vector<std::string> seq_tup, seq_names;
  kseq_t *seq;
  gzFile fp;
  int l;


  params.k_val_start = 5;
  params.k_val_end = 5;
  params.group_name = "sample";
  params.seq_filename = std::string("./") + params.group_name + std::string("/") + params.group_name + std::string(".pe.fastq.gz");
  params.hashes_filename = std::string("./reference/") + params.group_name + std::string(".hashes.bin");
  params.counts_filename = std::string("./reference/") + params.group_name + std::string(".counts.bin");
  params.indices_filename = std::string("./reference/") + params.group_name + std::string(".indices.bin");
  params.ids_filename = std::string("./reference/") + params.group_name + std::string(".ids.bin");

  params.num_threads = (params.k_val_end - params.k_val_start) / 2 + 1;
  params.is_ref = false;
  params.paired_end = true;

  KMerge *kmerge = new KMerge("lookup3", ".", "./reference", 0);

  params.kmerge = kmerge;

  KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
  tp.add_task(*task, &KMerge::BuilderTask::execute);
  tp.wait_for_all_tasks();


  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  uint seq_idx = 0;
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s), seq_name(seq->name.s);
    seq_tup.push_back(seq_str);
    seq_names.push_back(seq_name);
    if (seq_idx % 2 == 1) {
      std::string base_id = KMerge::get_seq_base_id(seq_names[0], seq_names[1]);
      for (uint k = params.k_val_start; k <= params.k_val_end; k+=2) {
	bool success = params.kmerge->hash_seq(seq_tup, k, hashed_counts, base_id, mtx);
	REQUIRE(success == true);
      }
      seq_tup.clear();
      seq_names.clear();
    }
    seq_idx++;
  }
  kseq_destroy(seq);
  gzclose(fp);


  delete kmerge;

  fs.open (params.hashes_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(hashes_in, ss);
  ss.str(std::string());

  fs.open (params.counts_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(counts_in, ss);
  ss.str(std::string());

  fs.open (params.indices_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(indices_in, ss);
  ss.str(std::string());

  for (uint i=0; i < indices_in.size() - 1; i++) {
    std::partial_sum(hashes_in.begin() + indices_in[i], hashes_in.begin() + indices_in[i+1], hashes_in.begin() + indices_in[i]);
  }


  fs.open (params.ids_filename.c_str(), ios::binary);
  decompressor.decompress(fs, ss);
  fs.close();
  dlib::deserialize(ids_in, ss);
  ss.str(std::string());

  REQUIRE(counts_in.size() == 118204);
  REQUIRE(hashes_in.size() == 118204);
  REQUIRE(indices_in.size() == 1001);
  REQUIRE(ids_in.size() == 1000);
  REQUIRE(indices_in.back() == 118204);
  uint total_hash_count=0;
  for (auto iter: hashed_counts) {
    total_hash_count += iter.second.size();
  }
  REQUIRE(total_hash_count == 118204);


  //ensure that hashes are in sorted order in each sample
  for (uint i=0; i < ids_in.size(); i++) {
    std::string seq_id = ids_in[i];
    uint last = 0;
    REQUIRE((indices_in[i+1]-indices_in[i]) == hashed_counts[seq_id].size());
    for (uint j=indices_in[i]; j < indices_in[i+1]; j++) {
      if (j != indices_in[i]) {
	REQUIRE(hashes_in[j] > last);
      }
      REQUIRE(hashed_counts[seq_id][hashes_in[j]] == counts_in[j]);
      last = hashes_in[j];
    }
  }

  if( remove( params.hashes_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

  if( remove( params.counts_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

  if( remove( params.indices_filename.c_str() ) != 0 )
    perror( "Error deleting file" );

  if( remove( params.ids_filename.c_str() ) != 0 )
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


  params.kmerge = new KMerge("lookup3", ".", "./reference", 0);
  bool success = params.kmerge->count_hashed_kmers(params, hashed_counts, true, true);
  REQUIRE(success == true);

  for(auto m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->second <= (uint) MAX_UINT_VAL);
  }
  hashed_counts.clear();
  delete params.kmerge;

  params.kmerge = new KMerge("spooky", ".", "./reference", 0);
  success = params.kmerge->count_hashed_kmers(params, hashed_counts, true, true);
  REQUIRE(success == true);

  for(auto m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    // make sure we are getting unsigned 32-bit integers back
    REQUIRE(m_iter->first >= 0);
    REQUIRE(m_iter->second <= (uint) MAX_UINT_VAL);
  }
  hashed_counts.clear();
  delete params.kmerge;

  params.kmerge = new KMerge("city", ".", "./reference", 0);
  success = params.kmerge->count_hashed_kmers(params, hashed_counts, true, true);
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
  std::string group("54095"); 
  std::string line;
  std::ifstream fs;
  std::map<std::string, std::string> taxonomy;


  
  KMerge *kmerge = new KMerge("lookup3", ".", "./reference", 0);

  kmerge->add_taxonomy(group);

  delete kmerge;

  fs.open (std::string("./reference/" + group + ".taxonomy.bin").c_str(), ios::binary);
  dlib::deserialize(taxonomy, fs);
  fs.close();

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

  if( remove(std::string("./reference/" + group + ".taxonomy.bin").c_str() ) != 0 )
    perror( "Error deleting file" );
  

}

