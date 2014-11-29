#include "kmerge.h"
#include "lookup3.c"
#include "MurmurHash3.h"
#include "city.h"
#include <dlib/serialize.h>
#include <fstream>
#include <unistd.h>

using namespace std;

KMerge::KMerge (const std::string& hash_func, const std::string& dir, const std::string& out_dir): dlog("kmerge") {
  this->dir = dir;
  this->out_dir = out_dir;
  this->dlog.set_level(dlib::LALL);

  if (hash_func == "lookup3") {
    this->hash_type = LOOKUP3;
  } else if (hash_func == "spooky") {
    this->hash_type = SPOOKY;
  } else if (hash_func == "murmur") {
    this->hash_type = MURMUR;
  } else if (hash_func == "city") {
    this->hash_type = CITY;
  } else {
    throw "Invalid hash function provided";
  }

}


KMerge::~KMerge() {
}

std::string KMerge::rev_comp(const std::string& input) {
  std::string output;

  for (std::string::const_reverse_iterator rit=input.rbegin(); rit!=input.rend(); rit++) {
    if (*rit == 'A') {
      output.append("T");
    } else if (*rit == 'T') {
      output.append("A");
    } else if (*rit == 'G') {
      output.append("C");
    } else if (*rit == 'C') {
      output.append("G");
    } else {
      throw "Invalid character in input";
    }
  }

  return output;
}

uint KMerge::hash_kmer(const std::string& kmer, const HashEnumType hash_type) {
  std::string rc_kmer = KMerge::rev_comp(kmer), combine;
  if (kmer < rc_kmer) {
    combine = kmer + rc_kmer;
  } else {
    combine = rc_kmer + kmer;
  }

  switch (hash_type) {
  case LOOKUP3:
    return(hashlittle(combine.c_str(), combine.length(), 0));
  case SPOOKY:
      return(SpookyHash::Hash32(combine.c_str(), combine.length(), 0));
  case MURMUR:
    uint out;
    MurmurHash3_x86_32  ( combine.c_str(), combine.length(), 0, &out );
    return out;
  case CITY:
    return(CityHash32(combine.c_str(), combine.length()));
  default:
    throw "Invalid hash function provided";
  }
}

uint KMerge::hash_kmer(const std::string& kmer) {
  return KMerge::hash_kmer(kmer, this->hash_type);
}


bool KMerge::count_hashed_kmers(param_struct& params,  btree::btree_map<uint, uint>& hashed_counts, bool split, bool print_status) {
  std::vector<std::tuple<uint, uint, uint, uint> > coords;
  uint pieces, piece_length;
  std::mutex mtx;
  kseq_t *seq;
  gzFile fp;
  std::vector<std::string> seqs;
  uint seq_id;
  int l;

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
	if (piece_length < k) {
	  piece_length = k;
	}
	if (!split) {
	  piece_length = str_len; // if hashing FASTQ sequences, don't split
	}
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
  KMerge::HashSeq func(params, seqs, hashed_counts, coords, mtx, print_status);
  dlib::parallel_for(params.num_threads, 0, coords.size(), func);
  coords.clear();
  kseq_destroy(seq);
  gzclose(fp);

  return true;
}

bool KMerge::hash_seq(std::string& seq, uint k, btree::btree_map<uint, uint>& hashed_counts, std::mutex& mtx) {
  std::vector<uint> hashes;
  uint hash;

  for (uint j = 0; j < seq.size() - k + 1; j++) {
    std::string kmer = seq.substr(j, k);
    std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
    if(kmer.find_first_not_of("ACGT") != std::string::npos) { 
      // skip kmers containing non-nucleotides            
      continue;
    }
    hash = this->hash_kmer(kmer);
    hashes.push_back(hash);
  }

  mtx.lock();
  for (auto it = hashes.begin(); it != hashes.end(); it++) hashed_counts[*it]++;
  mtx.unlock();

  hashes.clear();
  std::vector<uint>().swap(hashes);
  return true;
}


 
bool KMerge::add_dataset(const std::vector<uint>& data, const std::string& filename){
  cs compressor;
  std::stringstream ss;
  std::ofstream fs(filename.c_str(), ios::binary);

  dlib::serialize(data, ss);
  compressor.compress(ss, fs);
  fs.close();
  ss.str(std::string());

  return true;
}


bool KMerge::add_taxonomy(const std::string& group_name) {
  std::map<std::string, std::string> taxonomy;
  std::stringstream path_root, in_file_ss;
  std::string line;
  std::ofstream fs(std::string(this->out_dir + "/" + group_name + ".taxonomy.bin").c_str(), ios::binary);



  this->dlog << dlib::LINFO << "Adding taxonomy for " << group_name;

  path_root << "/" << group_name << "/taxonomy";

  in_file_ss << this->dir << path_root.str() << ".txt";
  ifstream in_file(in_file_ss.str().c_str());


  while (std::getline(in_file, line)) {
    std::istringstream tokenizer(line);
    std::string taxon, classification;
    uint i = 0;
    while (!tokenizer.eof()) {
      std::string token;
      getline(tokenizer, token, '\t');
      if (i == 0) {
	taxon = token;
      } else {
	classification = token;
      }
      i++;
    }
    taxonomy[taxon] = classification;
  }
  
  in_file.close();
  dlib::serialize(taxonomy, fs);

  this->dlog <<dlib::LINFO << "Finished adding taxonomy for " << group_name;

  return true;
}

 
void KMerge::BuilderTask::execute() {
  uint nz_count;
  btree::btree_map<uint, uint> hashed_counts;

  params.kmerge->dlog << dlib::LINFO << "Working on " << params.group_name;
  if(!(params.kmerge->count_hashed_kmers(params, hashed_counts, params.is_ref, params.is_ref))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to parse " << params.seq_filename;
    return;
  } else {
    params.kmerge->dlog << dlib::LINFO << "Finished parsing: " << params.seq_filename;  
  }

  nz_count = hashed_counts.size();

  params.kmerge->dlog << dlib::LINFO << "Hashes vector size: " << nz_count << " for " << params.group_name;


  params.kmerge->dlog <<dlib::LINFO << "Separating hashes and counts for " << params.group_name;
  std::vector<uint> hashes(nz_count), counts(nz_count);

  uint i = 0;
  uint last;
  for (auto m_iter = hashed_counts.begin(); m_iter != hashed_counts.end(); m_iter++) {
    if (i == 0) {
      hashes[i] = m_iter->first;
      last = m_iter->first;
    } else {
      hashes[i] = m_iter->first - last;
      last = m_iter->first;
    }
    counts[i] = m_iter->second;
    i++;
  }

  params.kmerge->dlog <<dlib::LINFO << "Finished separating hashes and counts for " << params.group_name;

  // remove all elements from map as they are no longer needed
  hashed_counts.clear();
  btree::btree_map<uint,uint>().swap(hashed_counts);

  params.kmerge->dlog <<dlib::LINFO << "Writing data for " << params.group_name;

  if(!params.kmerge->add_dataset(hashes, params.hashes_filename)) {                                                
    params.kmerge->dlog << dlib::LERROR << "Unable to add hashes for " << params.group_name;
    return;                                                                                                                                                                          
  }                                                                
  
  hashes.clear();
  std::vector<uint>().swap( hashes );


  if(!params.kmerge->add_dataset(counts, params.counts_filename)) {                                                     
    params.kmerge->dlog << dlib::LERROR << "Unable to add counts for " << params.group_name;
    return;                                                                                                                                                                          
  }

  counts.clear();
  std::vector<uint>().swap( counts );
    
  if(params.is_ref) {
    if(!(params.kmerge->add_taxonomy(params.group_name))) {
      params.kmerge->dlog << dlib::LERROR << "Unable to add classifications for " << params.group_name;
      return;
    }
  }

  params.kmerge->dlog << dlib::LINFO << "Done (" << params.group_name  << ")";

  return;

}

void KMerge::HashSeq::operator() (long i) const {
  std::tuple<uint, uint, uint, uint> coord = coords[i];
  uint seq_id = std::get<0>(coord);
  std::string seq = seqs[seq_id];
  uint k = std::get<1>(coord);
  uint start = std::get<2>(coord);
  uint end = std::get<3>(coord);
  std::string piece(seq.begin() + start, seq.begin() + end);

  if (print_status) params.kmerge->dlog << dlib::LINFO << "Processing " << params.group_name << "|" << seq_id << "|" << start << ":" << end;
  
  params.kmerge->hash_seq(piece, k, hashed_counts, mtx);
 
  if (print_status) params.kmerge->dlog << dlib::LINFO << "Finished processing " << params.group_name << "|" << seq_id << "|" << start << ":" << end;
}

