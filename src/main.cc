#include <dirent.h>
#include <sys/stat.h>
#include <string>
#include <fstream>
#include "kmerge.h"
#include <dlib/threads.h>
#include <dlib/cmd_line_parser.h>

using namespace std;

int main(int argc, char const ** argv) {
  dlib::command_line_parser parser;

  parser.add_option("o","Hash database file name", 1);
  parser.add_option("k", "Start and end k-mer values", 2);
  parser.add_option("d", "Location of sequences and taxonomy directories", 1);
  parser.add_option("i", "Location of single sequence file", 1);
  parser.add_option("t", "Max number of threads to use", 1);
  parser.add_option("f", "Hash function to use for k-mers", 1);
  parser.add_option("p", "Lock priority for database file name", 1);
  parser.add_option("g", "Max memory (in GB)", 1);
  parser.add_option("h","Display this help message.");


  //Parse command line.
  parser.parse(argc,argv);

  if (!(parser.parsed_line())) {
    std::cerr << "Line not parsed" << std::endl;
    return 0;
  }
  //Ensure options only defined once
  const char* one_time_opts[] = {"o", "k", "d", "i", "t", "f", "p", "g", "h"};
  parser.check_one_time_options(one_time_opts);

  //Check only one of d and i options given
  parser.check_incompatible_options("d", "i");

  parser.check_option_arg_range("k", 3, 31);
  parser.check_option_arg_range("t", 1, 80);

  // check if the -h option was given on the command line
  if (parser.option("h")) {
    // display all the command line options
    std::cout << "Usage: kmerge -o output_file -k k_start k_end (-d|-i|-t|-f)" << std::endl;
    parser.print_options(); 
    return 0;
  }

  // make sure two values input to k
  if (parser.option("k").number_of_arguments() != 2) {
    std::cerr << "Error in command line:\n   You must specify start and end k-mer values for parsing" << std::endl;
    std::cerr << "\nTry the -h option for more information." << std::endl;
    return 0;
  }
 
  // make sure one of the c or d options was given
  if (!parser.option("o")) {
    std::cerr << "Error in command line:\n   You must specify an database output filename" << std::endl;
    std::cerr << "\nTry the -h option for more information." << std::endl;
    return 0;
  }


  uint k_val_start = 0;
  uint k_val_end = 0;
  std::string db_filename, seq_dir("."), hash_func("lookup3"), in_file;
  uint num_threads = 1, max_threads = 0, parallel_for_threads = 1, priority = 1;
  double max_gb = 70;

  db_filename = parser.option("o").argument();
  k_val_start = atoi(parser.option("k").argument(0).c_str());
  k_val_end = atoi(parser.option("k").argument(1).c_str());


  if ((k_val_start % 2 == 0) || (k_val_end % 2 == 0)) {
    std::cerr << "Error in command line:\n   Start and end k-mer values must be odd" << std::endl;
    std::cerr << "\nTry the -h option for more information." << std::endl;
    return 0;
  }
  if (parser.option("d")) {
    seq_dir = parser.option("d").argument();
  }
  if (parser.option("i")) {
    in_file = parser.option("i").argument();
  } 
  if (parser.option("t")) {
    max_threads = atoi(parser.option("t").argument().c_str());
    uint num_ks = (k_val_end - k_val_start) / 2 + 1;
    if (max_threads <= 1) {
      parallel_for_threads = 1;
    } else if (max_threads <= num_ks) { //need 1 thread allocated for calling thread and 1 for memory polling
      num_threads = 2;
      parallel_for_threads = max_threads - 2;
    } else { //need 1 thread allocated for calling thread and 1 for memory polling
      num_threads = (max_threads / (num_ks + 1)) * 2; //adding one to account for polling thread 
      parallel_for_threads = num_ks;
    }
  }
  if (parser.option("f")) {
    std::string option = parser.option("f").argument();
    // check that values are valid or use default
    if ((option == "lookup3") || (option == "spooky") || (option == "murmur") || (option == "city")) {
      hash_func = option;
    }
  }
  if (parser.option("p")) {
    priority = atoi(parser.option("p").argument().c_str());
    if (priority < 0) {
      priority = 1;
    }
  }
  if (parser.option("g")) {
    max_gb = atoi(parser.option("g").argument().c_str());
  }

  stringstream seq_filename, file_loc, dataset_name;
  std::string group_name("");
  std::string hash_dataset_name("");
  std::string count_dataset_name("");
  KMerge *kmerge;

  std::vector<KMerge::BuilderTask*> task_ptrs;
  std::vector<ulib::chain_hash_map<uint, uint>*> map_ptrs;
  std::vector<std::mutex*> mtx_ptrs;
  std::vector<std::condition_variable*> cv_ptrs;
  dlib::thread_pool tp(num_threads);

  if (parser.option("i")) {
    kmerge = new KMerge(db_filename, hash_func, "", max_gb);
    param_struct params;
    params.kmerge = kmerge;
    params.db_filename = db_filename;
    params.lock_filename = params.db_filename + std::string(".lck");
    params.k_val_start = k_val_start;
    params.k_val_end = k_val_end;
    params.group_name = "sample";
    seq_filename.str("");
    params.seq_filename = in_file;
    dataset_name.str("");
    params.num_threads = parallel_for_threads;
    params.priority = priority;
    params.dump_mtx = new std::mutex;
    mtx_ptrs.push_back(params.dump_mtx);
    params.cv = new std::condition_variable;
    cv_ptrs.push_back(params.cv);
    params.ready = true;
    params.finished_hashing = false;
    params.hashed_counts = new ulib::chain_hash_map<uint, uint>(100000000);
    map_ptrs.push_back(params.hashed_counts);
    params.is_ref = false;
    params.dump_filename = "sample_hashes.bin";
    params.writing.resize(params.num_threads, 0);
    params.polling_done = true;
    KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
    if (max_threads > 1) {
      params.polling_done = false;
      tp.add_task(*task, &KMerge::BuilderTask::check_memory);
    }
    tp.add_task(*task, &KMerge::BuilderTask::execute);
    task_ptrs.push_back(task);
 
  } else {
    struct stat st;
    DIR *dirp;
    struct dirent *dp;
    dirp = opendir(seq_dir.c_str());
    kmerge = new KMerge(db_filename, hash_func, seq_dir, max_gb);
    while ((dp = readdir(dirp)) != NULL) {
      if(stat(dp->d_name, &st) == 0) {
	if (S_ISDIR(st.st_mode)) {
	  std::string s_org(dp->d_name);
	  if (s_org.compare(".") != 0 && s_org.compare("..") != 0) {
	    try {
	      param_struct params;
	      params.kmerge = kmerge;
	      params.db_filename = db_filename;
	      params.lock_filename = params.db_filename + std::string(".lck");
	      params.k_val_start = k_val_start;
	      params.k_val_end = k_val_end;
	      params.group_name = s_org;
	      seq_filename.str("");
	      seq_filename << seq_dir << "/" << s_org << "/" << s_org << ".fasta.gz";
	      params.seq_filename = seq_filename.str();
	      params.num_threads = parallel_for_threads;
	      params.priority = priority;
	      params.dump_mtx = new std::mutex;
	      mtx_ptrs.push_back(params.dump_mtx);
	      params.cv = new std::condition_variable;
	      cv_ptrs.push_back(params.cv);
	      params.ready = true;
	      params.finished_hashing = false;
	      params.hashed_counts = new ulib::chain_hash_map<uint, uint>(100000000);
	      map_ptrs.push_back(params.hashed_counts);
	      params.is_ref = true;
	      params.dump_filename = s_org + std::string("_hashes.bin");
	      params.writing.resize(params.num_threads, 0);
	      params.polling_done = true;
	      KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
	      if (max_threads > 1) {
		params.polling_done = false;
		tp.add_task(*task, &KMerge::BuilderTask::check_memory);
	      }
	      tp.add_task(*task, &KMerge::BuilderTask::execute);
 	      task_ptrs.push_back(task);
	    } catch( exception &e ) {
	      cout << e.what() << std::endl;
	      return 1;
	    }
	  }
	}
      }
    }

    (void)closedir(dirp);
  }

  tp.wait_for_all_tasks();
  for (vector<KMerge::BuilderTask*>::iterator iter=task_ptrs.begin(); iter!=task_ptrs.end(); iter++) {
    delete *iter;
  }

  for (std::vector<ulib::chain_hash_map<uint, uint>*>::iterator iter=map_ptrs.begin(); iter!=map_ptrs.end(); iter++) {
    delete *iter;
  }

  for (std::vector<std::mutex*>::iterator iter=mtx_ptrs.begin(); iter!=mtx_ptrs.end(); iter++) {
    delete *iter;
  }

  for (std::vector<std::condition_variable*>::iterator iter=cv_ptrs.begin(); iter!=cv_ptrs.end(); iter++) {
    delete *iter;
  }


  delete kmerge;
  

  return 0;
}
