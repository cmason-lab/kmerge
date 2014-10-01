#include "kmerge.h"
#include "lookup3.c"
#include "MurmurHash3.h"
#include "city.h"
#include <dlib/serialize.h>
#include <fstream>
#include <unistd.h>
#include "fastpfor/codecfactory.h"
#include "getRSS.c"

using namespace std;

pthread_mutex_t KMerge::db_mutex = PTHREAD_MUTEX_INITIALIZER;

KMerge::KMerge (const std::string& filename, const std::string& hash_func, const std::string& dir, double max_gb): dlog("kmerge") {
  this->filename = filename;
  this->dir = dir;
  this->dlog.set_level(dlib::LALL);
  if (max_gb < MIN_MEM) {
    this->dlog << dlib::LINFO << "Setting minimum memory to " << MIN_MEM << " GB";
    this->max_gb = MIN_MEM;
  } else {
    this->max_gb = max_gb;
  }

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

double KMerge::memory_used(void) {
  return getCurrentRSS() / BYTES_IN_GB;
}

void KMerge::dump_hashes(ulib::chain_hash_map<uint, uint>& o_map, std::string& filename) {
  std::ofstream ofs(filename.c_str(), ios::out | ios::app | ios::binary);
  for (auto it = o_map.begin(); it != o_map.end(); it++) {
    ofs << it.key() << "\t" << it.value() << std::endl;
  }
  o_map.clear();
  ofs.close();
}

void KMerge::load_hashes(btree::btree_map<uint, uint>& i_map, std::string& filename) {
  std::ifstream ifs(filename.c_str(), ios::in | ios::binary);
  uint hash, value;

  while (ifs >> hash >> value) {
    i_map[hash] += value;
  }

  ifs.close();
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
  params.finished_hashing = true; // stop memory polling
  return true;
}

std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > KMerge::compress(const std::vector<uint>& data) {
  FastPForLib::IntegerCODEC & codec =  * FastPForLib::CODECFactory::getFromName("simdfastpfor");
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > compressed(data.size()*1.1);
  size_t compressedsize = compressed.size();
  codec.encodeArray(data.data(), data.size(),
                    compressed.data(), compressedsize);
  compressed.resize(compressedsize);
  compressed.shrink_to_fit();
  return compressed;
}

  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > KMerge::uncompress(const std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> >& compressed, uint uncompressed_length) {
  FastPForLib::IntegerCODEC & codec =  * FastPForLib::CODECFactory::getFromName("simdfastpfor");
  std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > data(uncompressed_length);
  size_t recoveredsize = data.size();
  codec.decodeArray(compressed.data(),
                    compressed.size(), data.data(), recoveredsize);
  data.resize(recoveredsize);
  return data;

}

std::vector<uint> KMerge::compress2(const std::vector<uint>& data) {
  int num_bytes = data.size() * sizeof(uint);
  char* dest = new char[num_bytes];
  int comp_size = LZ4_compress((char*) data.data(), dest, num_bytes);
  std::vector<uint> compressed((uint*)dest, (uint*)dest + comp_size/sizeof(uint));
  delete [] dest;
  return compressed;
}

std::vector<uint> KMerge::uncompress2(const std::vector<uint>& compressed, uint uncompressed_length) {
  int num_bytes = uncompressed_length * sizeof(uint);
  char* dest = new char[num_bytes];
  int comp_size = LZ4_decompress_fast ((char*) compressed.data(), dest, num_bytes);
  std::vector<uint> data((uint*)dest, (uint*)dest + num_bytes/sizeof(uint));
  delete [] dest;
  return data;

}



bool KMerge::add_property(uint size, const std::string& key, unqlite *db) {

  this->dlog << dlib::LINFO << "Adding property " << key;

  int rc = unqlite_kv_store(db,key.c_str(),-1,&size,sizeof(uint));
  if (rc != UNQLITE_OK) {
    this->dlog << dlib::LERROR << "Unable to add property " << key << " - Error code: " << rc;
    return false;
  }

  this->dlog <<dlib::LINFO << "Finished adding property " << key;

  return true;

}
 


bool KMerge::add_dataset(const std::vector<uint>& data, const std::string& key, unqlite* db, bool compress){
  int rc;

  this->dlog << dlib::LINFO << "Adding data for " << key;

  if (compress == true) {
    this->dlog << dlib::LINFO << "Compressing data for " << key;
    std::vector<uint, FastPForLib::AlignedSTLAllocator<uint, BYTE_ALIGNED_SIZE> > compressed = KMerge::compress(data);
    this->dlog << dlib::LINFO << data.size() << "-->" << compressed.size();
    this->dlog << dlib::LINFO << "Finished compressing data for " << key;
    rc = unqlite_kv_store(db,key.c_str(),-1,&compressed[0],sizeof(uint)*compressed.size());
    this->add_property(compressed.size(), (key + "|size").c_str(), db);
  } else {
    rc = unqlite_kv_store(db,key.c_str(),-1,&data[0],sizeof(uint)*data.size());
    this->add_property(data.size(), (key + "|size").c_str(), db);
  }
  if (rc != UNQLITE_OK) {
    this->dlog << dlib::LERROR << "Unable to write data for " << key << " - Error code: " << rc;
    return false;
  }

  this->dlog <<dlib::LINFO << "Finished adding data for " << key;

  return true;

}



bool KMerge::add_taxonomy(const std::string& group_name, unqlite* db) {
  std::map<std::string, std::string> taxonomy;
  std::stringstream ss_out, path_root, in_file_ss;
  std::string line;



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
  dlib::serialize(taxonomy, ss_out);

  int rc = unqlite_kv_store(db,(group_name + std::string("|taxonomy")).c_str(),-1, ss_out.str().c_str(), ss_out.str().size());

  if (rc != UNQLITE_OK) {
    this->dlog << dlib::LERROR << "Unable to write taxonomy data for " << group_name << " - Error code: " << rc;
    return false;
  }
  this->add_property(ss_out.str().size(), (group_name + std::string("|taxonomy|size")).c_str(), db);


  this->dlog <<dlib::LINFO << "Finished adding taxonomy for " << group_name;

  return true;
}

void KMerge::poll_memory(param_struct& params) {
  while (params.finished_hashing == false) {
    if (KMerge::memory_used() > params.kmerge->max_gb) {
      std::unique_lock<std::mutex> lck(*(params.dump_mtx));
      params.ready = false;
      while (!(std::all_of(params.writing.begin(), params.writing.end(), [](char c){return c == 0;}))) sleep(1); // make sure no thread is writing to hashed_counts before proceeding
      KMerge::dump_hashes(*(params.hashed_counts), params.dump_filename);
      params.ready = true;
      params.cv->notify_all();
    }
    sleep(POLL_INTERVAL);
  }
}

void KMerge::build(param_struct& params) {
  stringstream file_name, file_loc;
  ulib::chain_hash_map<uint, uint> hashed_counts(KMerge::INIT_MAP_CAPACITY);

  params.kmerge->dlog << dlib::LINFO << "Working on " << params.group_name;
  if(!(params.kmerge->count_hashed_kmers(params, hashed_counts, false))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to parse " << params.seq_filename;
    return;
  } else {
    params.kmerge->dlog << dlib::LINFO << "Finished parsing: " << params.seq_filename;  
  }

  params.finished_hashing = true;

  uint nz_count = hashed_counts.size();

  params.kmerge->dlog << dlib::LINFO << "Hashes vector size: " << nz_count;

  params.kmerge->dlog <<dlib::LINFO << "Sorting hashes for " << params.group_name;


  std::vector<uint> hashes, counts;

  for (ulib::chain_hash_map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    hashes.push_back(iter.key());
  }
  std::sort(hashes.begin(), hashes.end());

  params.kmerge->dlog <<dlib::LINFO << "Finished sorting hashes for " << params.group_name;

  params.kmerge->dlog <<dlib::LINFO << "Separating hashes and counts for " << params.group_name;


  for (std::vector<uint>::const_iterator h_iter = hashes.begin(); h_iter != hashes.end(); ++h_iter) {
    counts.push_back(hashed_counts[*h_iter]);
  }

  params.kmerge->dlog <<dlib::LINFO << "Finished separating hashes and counts for " << params.group_name;
  
  // remove all elements from map as they are no longer needed
  hashed_counts.clear();


  int rc = unqlite_open(&(params.db), params.db_filename.c_str(), UNQLITE_OPEN_CREATE);
  if (rc != UNQLITE_OK) {
    params.kmerge->dlog << dlib::LERROR << params.group_name << " - Error code: " << rc;
    return;
  }

  uint partitions = ceil((double) hashes.size()/KMerge::PARTITION_SIZE);

  if(!params.kmerge->add_property(partitions, params.group_name + std::string("|parts"), params.db)) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add dataset size for " << params.group_name;
    unqlite_close(params.db);
    return;
  }

  for (uint i = 0; i < partitions; i++) {
    uint start_offset = i*KMerge::PARTITION_SIZE;
    uint end_offset = (i+1)*KMerge::PARTITION_SIZE;
    if (nz_count - i*KMerge::PARTITION_SIZE < KMerge::PARTITION_SIZE) end_offset = nz_count;
    std::vector<uint> sub_hashes(hashes.begin() + start_offset, hashes.begin() + end_offset);
    std::vector<uint> sub_counts(counts.begin() + start_offset, counts.begin() + end_offset);
    if(!params.kmerge->add_dataset(sub_hashes, params.group_name + std::string("|kmer_hash|") + std::to_string(i), params.db, false)) {
      params.kmerge->dlog << dlib::LERROR << "Unable to add hashes for " << params.group_name;
      unqlite_close(params.db);
      return;
    }
    sub_hashes.clear();
    std::vector<uint>().swap(sub_hashes);

    if(!params.kmerge->add_dataset(sub_counts, params.group_name + std::string("|count|") + std::to_string(i), params.db, true)) {
      params.kmerge->dlog << dlib::LERROR << "Unable to add counts for " << params.group_name;
      unqlite_close(params.db);
      return;
    }
    sub_counts.clear();
    std::vector<uint>().swap(sub_counts);
  }
  unqlite_close(params.db);

  params.kmerge->dlog << dlib::LINFO << hashes.size() << " k-mer hashes for " << params.group_name;

  params.kmerge->dlog << dlib::LINFO << "Done (" << params.group_name  << ")";

  std::vector<uint>().swap( hashes );
  std::vector<uint>().swap( counts );

  if( remove( params.dump_filename.c_str() ) != 0 )
    params.kmerge->dlog << dlib::LERROR << "Error deleting dump file";


  return;

}

void KMerge::BuilderTask::check_memory() {
  params.kmerge->poll_memory(params);
  params.polling_done = true; // allow execute to proceed
}

void KMerge::BuilderTask::execute() {
  stringstream file_name, file_loc;
  uint nz_count;

  params.kmerge->dlog << dlib::LINFO << "Working on " << params.group_name;
  if(!(params.kmerge->count_hashed_kmers(params, *(params.hashed_counts), params.is_ref))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to parse " << params.seq_filename;
    return;
  } else {
    params.kmerge->dlog << dlib::LINFO << "Finished parsing: " << params.seq_filename;  
  }

  while (!params.polling_done) sleep(5); 
  
  KMerge::dump_hashes(*(params.hashed_counts), params.dump_filename);

  btree::btree_map<uint, uint> tmp;
  
  params.kmerge->dlog <<dlib::LINFO << "Loading hashes for " << params.group_name;
  KMerge::load_hashes(tmp, params.dump_filename);
  params.kmerge->dlog <<dlib::LINFO << "Finished loading hashes for " << params.group_name;

  nz_count = tmp.size();
  params.kmerge->dlog << dlib::LINFO << "Hashes vector size: " << nz_count << " for " << params.group_name;


  params.kmerge->dlog <<dlib::LINFO << "Separating hashes and counts for " << params.group_name;
  std::vector<uint> hashes(nz_count), counts(nz_count);

  uint i = 0;
  for (btree::btree_map<uint, uint>::const_iterator m_iter = tmp.begin(); m_iter != tmp.end(); m_iter++) {
    hashes[i] = m_iter->first;
    counts[i] = m_iter->second;
    i++;
  }

  params.kmerge->dlog <<dlib::LINFO << "Finished separating hashes and counts for " << params.group_name;

  // remove all elements from map as they are no longer needed
  tmp.clear();
  btree::btree_map<uint,uint>().swap(tmp);


  while (mkdir(params.lock_filename.c_str(), 0644) == -1) sleep(params.priority); //process level lock

  
  params.kmerge->dlog << dlib::LINFO << params.group_name << " obtained db lock";

  int rc = unqlite_open(&(params.db), params.db_filename.c_str(), UNQLITE_OPEN_CREATE);
  if (rc != UNQLITE_OK) {
    params.kmerge->dlog << dlib::LERROR << params.group_name << " - Error code: " << rc;
    return;
  }

  params.kmerge->add_property(nz_count, (params.group_name + std::string("|size")).c_str(), params.db);

  uint partitions = ceil((double) hashes.size()/KMerge::PARTITION_SIZE);

  if(!params.kmerge->add_property(partitions, params.group_name + std::string("|parts"), params.db)) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add dataset size for " << params.group_name;
    unqlite_close(params.db);
    return;
  }


  for (uint i = 0; i < partitions; i++) {
    uint start_offset = i*KMerge::PARTITION_SIZE;
    uint end_offset = (i+1)*KMerge::PARTITION_SIZE;
    if (nz_count - i*KMerge::PARTITION_SIZE < KMerge::PARTITION_SIZE) end_offset = nz_count;
    std::vector<uint> sub_hashes(hashes.begin() + start_offset, hashes.begin() + end_offset);
    std::vector<uint> sub_counts(counts.begin() + start_offset, counts.begin() + end_offset);
    if(!params.kmerge->add_dataset(sub_hashes, params.group_name + std::string("|kmer_hash|") + std::to_string(i), params.db, false)) {
      params.kmerge->dlog << dlib::LERROR << "Unable to add hashes for " << params.group_name;
      unqlite_close(params.db);
      return;
    }
    sub_hashes.clear();
    std::vector<uint>().swap(sub_hashes);

    if(!params.kmerge->add_dataset(sub_counts, params.group_name + std::string("|count|") + std::to_string(i), params.db, true)) {
      params.kmerge->dlog << dlib::LERROR << "Unable to add counts for " << params.group_name;
      unqlite_close(params.db);
      return;
    }
    sub_counts.clear();
    std::vector<uint>().swap(sub_counts);

  }
  if(params.is_ref) {
    if(!(params.kmerge->add_taxonomy(params.group_name, params.db))) {
      params.kmerge->dlog << dlib::LERROR << "Unable to add classifications for " << params.group_name;
      unqlite_close(params.db);
      return;
    }
  }


  unqlite_close(params.db);

  rmdir(params.lock_filename.c_str());

  params.kmerge->dlog <<dlib::LINFO << params.group_name << " relinquished db lock";

  params.kmerge->dlog << dlib::LINFO << hashes.size() << " k-mer hashes for " << params.group_name;

  params.kmerge->dlog << dlib::LINFO << "Done (" << params.group_name  << ")";

  
  hashes.clear();
  std::vector<uint>().swap( hashes );
  counts.clear();
  std::vector<uint>().swap( counts );

  if( remove( params.dump_filename.c_str() ) != 0 )
    params.kmerge->dlog << dlib::LERROR << "Error deleting dump file";


  return;

}

void KMerge::CountAndHashSeq::operator() (long i) const {
  uint hash, counter=0;
  int l;
  kseq_t *seq;
  gzFile fp;
  std::unique_lock<std::mutex> lck(*(params.dump_mtx));


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
      while (!params.ready) params.cv->wait(lck);
      params.writing[i] = 1;
      hashed_counts[hash]++;
      params.writing[i] = 0;
      counter++;
    }
  }
  gzclose(fp);
  params.kmerge->dlog << dlib::LINFO << "Parsed " << counter << " hashes from sequence(s) in " << params.group_name << " (k = " << k << ")";
}

