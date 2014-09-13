#include "kmerge.h"
#include "lookup3.c"
#include "MurmurHash3.h"
#include "city.h"
#include <dlib/serialize.h>
#include <fstream>
#include <unistd.h>
#include "fastpfor/codecfactory.h"

using namespace std;



KMerge::KMerge (const std::string& filename, const std::string& hash_func, const std::string& dir): dlog("kmerge") {
  this->filename = filename;
  this->dir = dir;
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

bool KMerge::add_hash_and_count(std::vector<uint>& hashes, std::vector<uint>& counts, uint kmer_hash_val, uint kmer_count) {
  std::vector<uint>::iterator iter;
  uint pos = 0;        

  for (iter=hashes.begin(); iter!=hashes.end(); iter++) { 
    if (*iter==kmer_hash_val) {
      break;
    }
    pos++;
  }
  if (iter == hashes.end()) {
    hashes.push_back(kmer_hash_val);
    counts.push_back(kmer_count);
  } else {
    counts[pos] += kmer_count;
  }


  return true;
}

bool KMerge::add_hash_and_count(std::map<uint, uint>& hashed_counts, uint kmer_hash_val, uint kmer_count) {
  if (hashed_counts.find(kmer_hash_val) == hashed_counts.end()) {
    hashed_counts[kmer_hash_val] = kmer_count;
  } else {
    hashed_counts[kmer_hash_val] += kmer_count;
  }
  return true;
}

bool KMerge::add_hash(btree::btree_map<uint, uint>& hashed_counts, uint kmer_hash_val) {
  if (hashed_counts.find(kmer_hash_val) == hashed_counts.end()) {
    hashed_counts[kmer_hash_val] = 1; 
  } else {   
    hashed_counts[kmer_hash_val]++; 
  }
  return true;
}

bool KMerge::count_hashed_kmers(param_struct& params,  ulib::chain_hash_map<uint, uint>& hashed_counts, bool print_status) {
  KMerge::CountAndHashSeq func(params, hashed_counts, print_status);
  dlib::parallel_for(params.num_threads, (params.k_val_start + 1) / 2, (params.k_val_end + 1) / 2 + 1, func);

  return true;
}

std::vector<uint> KMerge::compress(const std::vector<uint>& data) {
  FastPForLib::IntegerCODEC & codec =  * FastPForLib::CODECFactory::getFromName("simdfastpfor");
  std::vector<uint32_t> compressed(data.size()*1.1);
  size_t compressedsize = compressed.size();
  codec.encodeArray(data.data(), data.size(),
                    compressed.data(), compressedsize);
  compressed.resize(compressedsize);
  compressed.shrink_to_fit();
  return compressed;
}

std::vector<uint> KMerge::uncompress(const std::vector<uint>& compressed_output, uint uncompressed_length) {
  FastPForLib::IntegerCODEC & codec =  * FastPForLib::CODECFactory::getFromName("simdfastpfor");
  std::vector<uint32_t> data(uncompressed_length);
  size_t recoveredsize = data.size();

  codec.decodeArray(compressed_output.data(),
                    compressed_output.size(), data.data(), recoveredsize);
  data.resize(recoveredsize);

  return data;

}



bool KMerge::add_dataset_size(uint size, const std::string& key, leveldb::DB* db) {
  std::stringstream ss_out;

  ss_out << size;

  leveldb::Status s = db->Put(leveldb::WriteOptions(), key, ss_out.str());
  if (!s.ok()) {
    this->dlog << dlib::LERROR << "Unable to write data for " << key;
    return false;
  }
  return true;

}
 



bool KMerge::add_dataset(const std::vector<uint>& data, const std::string& key, leveldb::DB* db){
  std::stringstream ss_out;

  std::vector<uint> compressed = KMerge::compress(data);

  dlib::serialize(compressed, ss_out);

  leveldb::Status s = db->Put(leveldb::WriteOptions(), key, ss_out.str());
  if (!s.ok()) {
    this->dlog << dlib::LERROR << "Unable to write data for " << key;
    return false;
  }

  return true;

}



bool KMerge::add_taxonomy(const std::string& group_name, leveldb::DB* db) {
  std::map<std::string, std::string> taxonomy;
  std::stringstream ss_out, path_root, in_file_ss;
  std::string line;

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
  dlib::serialize(taxonomy, ss_out);

  leveldb::Status s = db->Put(leveldb::WriteOptions(), group_name + std::string("|taxonomy"), ss_out.str());
  if (!s.ok()) {
    this->dlog << dlib::LERROR << "Unable to write taxonomy data for " << group_name + std::string("|taxonomy");
    return false;
  }
  return true;
}



void KMerge::build(param_struct& params) {
  stringstream file_name, file_loc;
  ulib::chain_hash_map<uint, uint> hashed_counts(KMerge::INIT_MAP_CAPACITY);
  leveldb::Options options;
  options.create_if_missing = true;

  params.kmerge->dlog << dlib::LINFO << "Working on " << params.group_name;
  if(!(params.kmerge->count_hashed_kmers(params, hashed_counts, false))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to parse " << params.seq_filename;
    return;
  } else {
    params.kmerge->dlog << dlib::LINFO << "Finished parsing: " << params.seq_filename;  
  }

  params.kmerge->dlog << dlib::LINFO << "Hashes vector size: " << hashed_counts.size();

  btree::btree_map<uint, uint> sorted_hashed_counts;
  for (ulib::chain_hash_map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    sorted_hashed_counts[iter.key()] = iter.value();
  }

  // remove all elements from map as they are no longer needed
  hashed_counts.clear();

  std::vector<uint> hashes, counts;

  for (btree::btree_map<uint, uint>::iterator b_iter = sorted_hashed_counts.begin(); b_iter != sorted_hashed_counts.end(); ++b_iter) {
    hashes.push_back(b_iter->first);
    counts.push_back(b_iter->second);
  }

  
  // remove all elements from map as they are no longer needed
  sorted_hashed_counts.clear();
  btree::btree_map<uint, uint>().swap(sorted_hashed_counts);

  
  leveldb::Status s = leveldb::DB::Open(options, params.db_filename, &(params.db));
  if (!s.ok()) {
    params.kmerge->dlog << dlib::LERROR << s.ToString();
    return;
  }

  if(!params.kmerge->add_dataset_size(hashes.size(), params.group_name + std::string("|size"), params.db)) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add dataset size for " << params.group_name;
    return;
  }
  if(!params.kmerge->add_dataset(hashes, params.group_name + std::string("|kmer_hash"), params.db)) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add hashes for " << params.group_name;
    return;
  }

  if(!params.kmerge->add_dataset(counts, params.group_name + std::string("|count"), params.db)) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add counts for " << params.group_name;
    return;
  }

  delete params.db;

  params.kmerge->dlog << dlib::LINFO << hashes.size() << " k-mer hashes for " << params.group_name;

  params.kmerge->dlog << dlib::LINFO << "Done (" << params.group_name  << ")";

  std::vector<uint>().swap( hashes );
  std::vector<uint>().swap( counts );

  return;

}



void KMerge::BuilderTask::execute() {
  stringstream file_name, file_loc;
  ulib::chain_hash_map<uint, uint> hashed_counts(KMerge::INIT_MAP_CAPACITY);
  std::vector<uint> hashes, counts;
  leveldb::Options options;
  options.create_if_missing = true;

  params.kmerge->dlog << dlib::LINFO << "Working on " << params.group_name;
  if(!(params.kmerge->count_hashed_kmers(params, hashed_counts, true))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to parse " << params.seq_filename;
    return;
  } else {
    params.kmerge->dlog << dlib::LINFO << "Finished parsing: " << params.seq_filename;  
  }

  params.kmerge->dlog << dlib::LINFO << "Hashes vector size: " << hashed_counts.size() << " for " << params.group_name;

  // put hashes in sorted order via btree_map
  btree::btree_map<uint, uint> sorted_hashed_counts;
  for (ulib::chain_hash_map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    sorted_hashed_counts[iter.key()] = iter.value();
  }

  // remove all elements from map as they are no longer needed
  hashed_counts.clear();

  for (btree::btree_map<uint, uint>::iterator b_iter = sorted_hashed_counts.begin(); b_iter != sorted_hashed_counts.end(); ++b_iter) {
    hashes.push_back(b_iter->first);
    counts.push_back(b_iter->second);
  }
  
  // remove all elements from map as they are no longer needed
  sorted_hashed_counts.clear();
  btree::btree_map<uint, uint>().swap(sorted_hashed_counts);

  while (mkdir(params.lock_filename.c_str(), 0644) == -1) sleep(params.priority); //process level lock

  leveldb::Status s = leveldb::DB::Open(options, params.db_filename, &(params.db));
  if (!s.ok()) {
    params.kmerge->dlog << dlib::LERROR << s.ToString();
    return;
  }
  

  if(!params.kmerge->add_dataset_size(hashes.size(), params.group_name + std::string("|size"), params.db)){
    params.kmerge->dlog << dlib::LERROR << "Unable to add dataset size for " << params.group_name;
    return;
  }

  if(!params.kmerge->add_dataset(hashes, params.group_name + std::string("|kmer_hash"), params.db)) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add hashes for " << params.group_name;
    return;
  }

  if(!params.kmerge->add_dataset(counts, params.group_name + std::string("|count"), params.db)) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add counts for " << params.group_name;
    return;
  }

  if(!(params.kmerge->add_taxonomy(params.group_name, params.db))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add classifications for " << params.group_name;
    return;
  }

  delete params.db;
  
  rmdir(params.lock_filename.c_str());

  params.kmerge->dlog << dlib::LINFO << hashes.size() << " k-mer hashes for " << params.group_name;

  params.kmerge->dlog << dlib::LINFO << "Done (" << params.group_name  << ")";

  
  hashes.clear();
  std::vector<uint>().swap( hashes );
  counts.clear();
  std::vector<uint>().swap( counts );

  return;

}

void KMerge::CountAndHashSeq::operator() (long i) const {
  uint hash, counter=0;
  int l;
  kseq_t *seq;
  gzFile fp;


  uint k = 2*i - 1; // make sure k is converted to odd value from input index
  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    if (print_status == true) {
      params.kmerge->dlog << dlib::LINFO << "Processing " << params.group_name << " sequence " << seq->name.s << " for k = " << k;
    }
    std::string seq_str(seq->seq.s);
    if(seq_str.size() < k) return;
    for (uint j = 0; j < seq_str.size() - k + 1; j++) {
      std::string kmer = seq_str.substr(j, k);
      std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
      if(kmer.find_first_not_of("ACGT") != std::string::npos) { 
	// skip kmers containing non-nucleotides            
	continue;
      }
      hash = params.kmerge->hash_kmer(kmer);
      hashed_counts[hash]++;
      counter++;
    }
  }
  params.kmerge->dlog << dlib::LINFO << "Parsed " << counter << " hashes from sequence(s) in " << params.group_name << " (k = " << k << ")";
}

