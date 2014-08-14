#include "kmerge.h"
#include "lookup3.c"
#include "MurmurHash3.h"
#include "city.h"
#include <dlib/serialize.h>
#include <fstream>
#include <unistd.h>

using namespace std;

pthread_mutex_t KMerge::mutex = PTHREAD_MUTEX_INITIALIZER;



KMerge::KMerge (const std::string& filename, const std::string& hash_func, const std::string& dir): dlog("kmerge") {
  this->filename = filename;
  //this->hdf5_file = new HDF5(this->filename, false);
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
  //delete this->hdf5_file;
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

bool KMerge::count_hashed_kmers_fasta(param_struct& params, btree::btree_map<uint, uint>& hashed_counts) {
  int l;
  kseq_t *seq;
  gzFile fp;
  pthread_mutex_t mutex;

  pthread_mutex_init(&mutex, NULL);
  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string seq_str(seq->seq.s);
    KMerge::CountAndHashSeqFasta func(params, hashed_counts, seq_str, mutex);
    dlib::parallel_for(params.num_threads, (params.k_val_start + 1) / 2, (params.k_val_end + 1) / 2 + 1, func);
  }
  pthread_mutex_destroy(&mutex);
  kseq_destroy(seq);
  gzclose(fp);
  return true;
}


bool KMerge::count_hashed_kmers_fastq(param_struct& params, btree::btree_map<uint, uint>& hashed_counts) {
  pthread_mutex_t mutex;

  pthread_mutex_init(&mutex, NULL);
  KMerge::CountAndHashSeqFastq func(params, hashed_counts, mutex);
  dlib::parallel_for(params.num_threads, (params.k_val_start + 1) / 2, (params.k_val_end + 1) / 2 + 1, func);
  pthread_mutex_destroy(&mutex);

  return true;
}

bool KMerge::add_dataset(const uint data_size, const uint* data, param_struct& params) {
  H5File* file;
  H5::DataSet* dataset;
  uint num_cols = 1, rank;
  char *version, *date;
  int r, i;
  const std::string ds_name("counts");
  const hsize_t data_dims[2] = {data_size, num_cols};
  
  /* Register the filter with the library */
  r = register_blosc(&version, &date);

  ifstream ifile(this->filename);  
  if (ifile) { 
    file = new H5File( this->filename, H5F_ACC_RDWR );
    dataset = new H5::DataSet( file->openDataSet( ds_name ) );

    H5::DataSpace dataspace = dataset->getSpace();
    rank = dataspace.getSimpleExtentNdims();

    hsize_t dims_out[2];
    int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
    hsize_t offset[2];
    offset[0] = 0;
    offset[1] = dims_out[1];
    params.compound_name = params.group_name + std::string("|") + std::to_string(dims_out[1]);
    dims_out[1]++; // extend dataset by 1                                                                                                                                                                                                      
    dataset->extend( dims_out );

    DataSpace fspace = dataset->getSpace();
    fspace.selectHyperslab( H5S_SELECT_SET, data_dims, offset );
    DataSpace mspace( rank, data_dims );
    dataset->write( &data[0], H5::PredType::NATIVE_UINT, mspace, fspace );
    dataset->close();
    delete dataset;
  } else {
    file = new H5File( this->filename, H5F_ACC_TRUNC );
    unsigned int cd_values[7];
    /*
     * Create property list for a dataset and set up fill values. 
     */
    rank = 2;
    uint fillvalue = 0;   /* Fill value for the dataset */
    H5::DSetCreatPropList plist;
    hsize_t max_dims[2] = {data_size, H5S_UNLIMITED};
    hsize_t chunk_dims[2] = {KMerge::CHUNK_ROW_SIZE, num_cols};

    cd_values[4] = 1;       /* compression level */
    cd_values[5] = 1;       /* 0: shuffle not active, 1: shuffle active */
    cd_values[6] = BLOSC_LZ4HC; /* the actual compressor to use */

    plist.setChunk(rank, chunk_dims);
    plist.setFillValue(H5::PredType::NATIVE_UINT, &fillvalue);
    plist.setFilter(FILTER_BLOSC, H5Z_FLAG_OPTIONAL, 7, cd_values);
    plist.setShuffle();
    plist.setFletcher32();

    /* 
     * Create dataspace for the dataset 
     */
    H5::DataSpace fspace( rank, data_dims, max_dims);

    /*
     * Create dataset and write it into the file.
     */
    H5::DataSpace mspace( rank, data_dims, max_dims );
    dataset = new H5::DataSet( file->createDataSet(ds_name, H5::PredType::NATIVE_UINT, fspace, plist));
    dataset->write( &data[0], H5::PredType::NATIVE_UINT, mspace, fspace );
    dataset->close();
    delete dataset;
    params.compound_name = params.group_name + std::string("|0");
  }
  
  file->close();
  delete file;
  
  return true;
}

bool KMerge::add_dataset(const std::string dataset_path, const uint data_size, const uint* data, HDF5* h5_file=NULL) {
  std::vector<uint64_t> dims;

  HDF5 *hdf5_file = new HDF5(this->filename, false);

  dims.push_back(data_size);
  if (!(hdf5_file->createDataset(dataset_path, dims, FQ::FQT_UINT))) {
    this->dlog << dlib::LERROR << "Unable to create dataset " << dataset_path;
    return false;
  }
  if (!(hdf5_file->setData(dataset_path, data))) {
    this->dlog << dlib::LERROR << "Unable to set data for dataset" << dataset_path;
    return false;
  }

  delete hdf5_file;

  return true;
}


bool KMerge::add_taxonomy(const std::string& group_name, const std::string& compound_group_name) {
  std::stringstream path, in_file_ss, path_root, error;
  std::vector<std::string> lines;
  std::vector<uint64_t> dims;
  std::string line;

  HDF5* hdf5_file = new HDF5(this->filename, false);

  path_root << group_name << "/taxonomy";

  in_file_ss << this->dir << path_root.str() << ".txt";
  ifstream in_file(in_file_ss.str().c_str());

  dims.push_back(0);

  while (std::getline(in_file, line)) {
    std::istringstream tokenizer(line);
    path.str("");
    path << compound_group_name << "/taxonomy";
    while (!tokenizer.eof()) {
      std::string token;
      getline(tokenizer, token, '\t');
      path << "/" << token;
    }
    if (!(hdf5_file->createDataset(path.str(), dims, FQ::FQT_BYTE))) {
      this->dlog << dlib::LERROR << "Unable to add classification:" << path.str();
      return false;
    }
  }
  
  in_file.close();
  delete hdf5_file;

  return true;
}

void KMerge::build(param_struct& params) {
  stringstream file_name, file_loc;
  btree::btree_map<uint, uint> hashed_counts;

  params.kmerge->dlog << dlib::LINFO << "Working on " << params.group_name;
  if(!(params.kmerge->count_hashed_kmers_fastq(params, hashed_counts))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to parse " << params.seq_filename;
    return;
  } else {
    params.kmerge->dlog << dlib::LINFO << "Finished parsing: " << params.seq_filename;
    params.kmerge->dlog << dlib::LINFO << "Hashes vector size: " << hashed_counts.size();  
  }

  uint hash_count = hashed_counts.size();

  std::vector<uint> hashes, counts;

  for (btree::btree_map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    hashes.push_back(iter->first);
    counts.push_back(iter->second);
  }

  
  // remove all elements from map as they are no longer needed
  hashed_counts.clear();
  btree::btree_map<uint, uint>().swap( hashed_counts );


  if(!(params.kmerge->add_dataset(params.hash_dataset_name, hashes.size(), &hashes[0], NULL))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add hashes for " << params.group_name;
    pthread_mutex_unlock( &KMerge::mutex );
    return;
  }

  if(!(params.kmerge->add_dataset(params.counts_dataset_name, counts.size(), &counts[0], NULL))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add counts for " << params.group_name;
    pthread_mutex_unlock( &KMerge::mutex );
    return;
  }


  pthread_mutex_unlock( &KMerge::mutex );

  params.kmerge->dlog << dlib::LINFO << hash_count << " k-mer hashes for " << params.group_name;

  params.kmerge->dlog << dlib::LINFO << "Done (" << params.group_name  << ")";

  std::vector<uint>().swap( hashes );
  std::vector<uint>().swap( counts );

  return;

}



void KMerge::BuilderTask::execute() {
  stringstream file_name, file_loc;
  btree::btree_map<uint, uint> hashed_counts;

  params.kmerge->dlog << dlib::LINFO << "Working on " << params.group_name;
  if(!(params.kmerge->count_hashed_kmers_fasta(params, hashed_counts))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to parse " << params.seq_filename;
    return;
  } else {
    params.kmerge->dlog << dlib::LINFO << "Finished parsing: " << params.seq_filename;
    params.kmerge->dlog << dlib::LINFO << "Hashes vector size: " << hashed_counts.size();  
  }

  uint hash_count = hashed_counts.size();

  /*std::vector<uint> hashes, counts;

  for (btree::btree_map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {
    hashes.push_back(iter->first);
    counts.push_back(iter->second);
  }
  */
  
  std::vector<uint> counts(MAX_UINT_VAL);

  for (btree::btree_map<uint, uint>::iterator iter = hashed_counts.begin(); iter != hashed_counts.end(); ++iter) {                                                                                      
    counts[iter->first] = iter->second;  
  }                                                                                                                                                                                                      

  // remove all elements from map as they are no longer needed
  hashed_counts.clear();
  btree::btree_map<uint, uint>().swap( hashed_counts );
  /*
  ofstream out_hashes_file(params.tmp_hashes_filename.c_str(), ios::out | ios::binary), 
    out_counts_file(params.tmp_counts_filename.c_str(), ios::out | ios::binary);

  try {
    params.kmerge->dlog << dlib::LINFO << "Serializing hashes for  " << params.group_name; 
    dlib::serialize(hashes, out_hashes_file);
    params.kmerge->dlog << dlib::LINFO << "Finished serializing hashes for  " << params.group_name;
  } catch (dlib::serialization_error& e) {
    params.kmerge->dlog << dlib::LERROR << "Unable to serialize hashes (" << params.group_name << ")";
    return;
  }

  std::vector<uint>().swap( hashes );
  out_hashes_file.close();


  try {
    params.kmerge->dlog << dlib::LINFO << "Serializing counts for  " << params.group_name;
    dlib::serialize(counts, out_counts_file);
    params.kmerge->dlog << dlib::LINFO << "Finished serializing counts for  " << params.group_name;
  } catch (dlib::serialization_error& e) {
    params.kmerge->dlog << dlib::LERROR << "Unable to serialize counts (" << params.group_name << ")";
    return;
  }

  std::vector<uint>().swap( counts );
  out_counts_file.close();

  */
  //pthread_mutex_lock( &KMerge::mutex );
  /*
  ifstream in_hashes_file(params.tmp_hashes_filename.c_str(), ios::in | ios::binary), 
    in_counts_file(params.tmp_counts_filename.c_str(), ios::in | ios::binary);
  
  try {
    params.kmerge->dlog << dlib::LINFO << "De-serializing hashes for  " << params.group_name;
    dlib::deserialize(hashes, in_hashes_file);
    params.kmerge->dlog << dlib::LINFO << "Finished de-serializing hashes for  " << params.group_name;
  } catch (dlib::serialization_error& e) {
    params.kmerge->dlog << dlib::LERROR << "Unable to deserialize hashes (" << params.group_name << ")";
    pthread_mutex_unlock( &KMerge::mutex );
    return;
  }

  in_hashes_file.close();

  if(!(params.kmerge->add_dataset(params.hash_dataset_name, hashes.size(), &hashes[0], NULL))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add hashes for " << params.group_name;
    pthread_mutex_unlock( &KMerge::mutex );
    return;
  }
  std::vector<uint>().swap( hashes );
  try {
    params.kmerge->dlog << dlib::LINFO << "De-serializing counts for  " << params.group_name;
    dlib::deserialize(counts, in_counts_file);
    params.kmerge->dlog << dlib::LINFO << "Finished de-serializing counts for  " << params.group_name;
  } catch (dlib::serialization_error& e) {
    params.kmerge->dlog << dlib::LERROR << "Unable to deserialize counts (" << params.group_name << ")";
    pthread_mutex_unlock( &KMerge::mutex );
    return;
  }

  in_counts_file.close();
  */
  /*
  if(!(params.kmerge->add_dataset(params.counts_dataset_name, counts.size(), &counts[0], NULL))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add counts for " << params.group_name;
    pthread_mutex_unlock( &KMerge::mutex );
    return;
  }
  */
  
  while (mkdir(params.lock_filename.c_str(), 0644) == -1) sleep(params.priority);

  if(!(params.kmerge->add_dataset(counts.size(), &counts[0], params))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add counts for " << params.group_name;
    //pthread_mutex_unlock( &KMerge::mutex );
    return;
  }

  if(!(params.kmerge->add_taxonomy(params.group_name, params.compound_name))) {
    params.kmerge->dlog << dlib::LERROR << "Unable to add classifications for " << params.group_name;
    //pthread_mutex_unlock( &KMerge::mutex );
    return;
  }
  rmdir(params.lock_filename.c_str());
  
  counts.clear();
  std::vector<uint>().swap( counts );
  //pthread_mutex_unlock( &KMerge::mutex );

  params.kmerge->dlog << dlib::LINFO << hash_count << " k-mer hashes for " << params.group_name;

  params.kmerge->dlog << dlib::LINFO << "Done (" << params.group_name  << ")";

  /*
  remove(params.tmp_hashes_filename.c_str());
  remove(params.tmp_counts_filename.c_str());
  */
  return;

}

template<class InputIterT1, class InputIterT2, class OutputIterT, class Comparator, class Func>
OutputIterT merge_apply(
			InputIterT1 first1, InputIterT1 last1,
			InputIterT2 first2, InputIterT2 last2,
			OutputIterT result, Comparator comp, Func func) {
  while (true)
    {
      if (first1 == last1) return std::copy(first2, last2, result);
      if (first2 == last2) return std::copy(first1, last1, result);

      if (comp(*first1, *first2) < 0) {
	*result = *first1;
	++first1;
      } else if (comp(*first1, *first2) > 0) {
	*result = *first2;
	++first2;
      } else { 
	*result = func(*first1, *first2);
	++first1;
	++first2; 
      }   
      ++result;
    }
}

template<class T>
int compare_first(T a, T b) {
  return a.first - b.first;
}

template<class T>
T sum_pairs(T a, T b) {
  return std::make_pair(a.first, a.second + b.second);
}

void KMerge::CountAndHashSeqFasta::operator() (long i) const {
  btree::btree_map<uint, uint> local_map, result;
  uint hash;

  uint k = 2*i - 1; // make sure k is converted to odd value from input index 
  if(seq.size() < k) return;
  for (uint j = 0; j < seq.size() - k + 1; j++) {
    std::string kmer = seq.substr(j, k);
    std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
    if(kmer.find_first_not_of("ACGT") != std::string::npos) { // skip kmers containing non-nucleotides            
      continue;
    }
    hash = params.kmerge->hash_kmer(kmer);
    local_map[hash]++;
  }
  params.kmerge->dlog << dlib::LINFO << "Parsed " << local_map.size() << " hashes from sequence in " << params.group_name << " (k = " << k << ")";
  pthread_mutex_lock( &m );
  merge_apply(hashed_counts.begin(), hashed_counts.end(), local_map.begin(), local_map.end(), std::inserter(result, result.begin()), compare_first<pair<uint, uint> >, sum_pairs<pair<uint, uint> >);
  local_map.clear();
  btree::btree_map<uint, uint>().swap(local_map);
  result.swap(hashed_counts);
  result.clear();
  btree::btree_map<uint, uint>().swap(result);
  pthread_mutex_unlock( &m );
}

void KMerge::CountAndHashSeqFastq::operator() (long i) const {
  btree::btree_map<uint, uint> local_map, result;
  int l;
  uint hash;
  kseq_t *seq;
  gzFile fp;
  uint k = 2*i - 1; // make sure k is converted to odd value from input index 

  fp = gzopen(params.seq_filename.c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    if (seq->seq.l < k) return;
    std::string seq_str(seq->seq.s);
    for (int i = 0; i < seq->seq.l - k + 1; i++) {
      std::string kmer = seq_str.substr(i, k);
      std::transform(kmer.begin(), kmer.end(), kmer.begin(), ::toupper);
      if(kmer.find_first_not_of("ACGT") != std::string::npos) { // skip kmers containing non-nucleotides
        continue;
      }
      hash = params.kmerge->hash_kmer(kmer);
      local_map[hash]++;
    }
  }
  kseq_destroy(seq);
  gzclose(fp);


  params.kmerge->dlog << dlib::LINFO << "Parsed " << local_map.size() << " hashes from sequences in " << params.group_name << " (k = " << k << ")";

  pthread_mutex_lock( &m );
  merge_apply(hashed_counts.begin(), hashed_counts.end(), local_map.begin(), local_map.end(), std::inserter(result, result.begin()), compare_first<pair<uint, uint> >, sum_pairs<pair<uint, uint> >);
  local_map.clear();
  btree::btree_map<uint, uint>().swap(local_map);
  result.swap(hashed_counts);
  result.clear();
  btree::btree_map<uint, uint>().swap(result);
  pthread_mutex_unlock( &m );
}
