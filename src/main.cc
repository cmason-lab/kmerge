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

  parser.add_option("o","HDF5 file name", 1);
  parser.add_option("k", "Start and end k-mer values", 2);
  parser.add_option("d", "Location of sequences and taxonomy directories", 1);
  parser.add_option("i", "Location of single sequence file", 1);
  parser.add_option("t", "Max number of threads to use", 1);
  parser.add_option("f", "Hash function to use for k-mers", 1);
  parser.add_option("h","Display this help message.");


  //Parse command line.
  parser.parse(argc,argv);

  if (!(parser.parsed_line())) {
    std::cerr << "Line not parsed" << std::endl;
    return 0;
  }
  //Ensure options only defined once
  const char* one_time_opts[] = {"o", "k", "d", "i", "t", "f", "h"};
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
    std::cerr << "Error in command line:\n   You must specify an hdf5 output file" << std::endl;
    std::cerr << "\nTry the -h option for more information." << std::endl;
    return 0;
  }


  uint k_val_start = 0;
  uint k_val_end = 0;
  std::string hdf5_filename, seq_dir("."), hash_func("lookup3"), in_file;
  uint num_threads = 1, max_threads = 0, parallel_for_threads = 1;

  hdf5_filename = parser.option("o").argument();
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
    if (max_threads == 1) {
      parallel_for_threads = 1;
    } else if (max_threads <= num_ks) { //need 1 thread allocated for calling thread 
      parallel_for_threads = max_threads - 1;
    } else {
      num_threads = max_threads / num_ks;
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
  
  stringstream seq_filename, file_loc, dataset_name, tmp_hashes_filename, tmp_counts_filename;
  vector<uint> hashes;
  vector<uint> counts;
  std::string group_name("");
  std::string hash_dataset_name("");
  std::string count_dataset_name("");

  KMerge *kmerge;

  if (parser.option("i")) {
    kmerge = new KMerge(hdf5_filename, hash_func, "");
    param_struct params;
    params.kmerge = kmerge;
    params.hdf5_filename = "sample.h5";
    params.k_val_start = k_val_start;
    params.k_val_end = k_val_end;
    params.group_name = "/sample";
    seq_filename.str("");
    params.seq_filename = in_file;
    dataset_name.str("");
    params.hash_dataset_name = "/sample/kmer_hash";
    params.counts_dataset_name = "/sample/count";
    params.tmp_hashes_filename = ""; 
    params.tmp_counts_filename = "";
    params.num_threads = parallel_for_threads;
    kmerge->build(params);
  } else {
    struct stat st;
    DIR *dirp;
    struct dirent *dp;
    vector<KMerge::BuilderTask*> task_ptrs;
    dlib::thread_pool tp(num_threads);
    dirp = opendir(seq_dir.c_str());
    kmerge = new KMerge(hdf5_filename, hash_func, seq_dir);
    while ((dp = readdir(dirp)) != NULL) {
      if(stat(dp->d_name, &st) == 0) {
	if (S_ISDIR(st.st_mode)) {
	  std::string s_org(dp->d_name);
	  if (s_org.compare(".") != 0 && s_org.compare("..") != 0) {
	    try {
	      param_struct params;
	      params.kmerge = kmerge;
	      params.hdf5_filename = hdf5_filename;
	      params.k_val_start = k_val_start;
	      params.k_val_end = k_val_end;
	      params.group_name = std::string("/") + s_org;
	      seq_filename.str("");
	      seq_filename << seq_dir << "/" << s_org << "/" << s_org << ".fasta.gz";
	      params.seq_filename = seq_filename.str();
	      dataset_name.str("");
	      dataset_name << "/" << s_org << "/" << "kmer_hash";
	      params.hash_dataset_name = dataset_name.str();
	      dataset_name.str("");
	      dataset_name << "/" << s_org << "/" << "count";
	      params.counts_dataset_name = dataset_name.str();
	      tmp_hashes_filename.str("");
	      tmp_hashes_filename << seq_dir << "/" << s_org << "/" << "hashes.bin";
	      params.tmp_hashes_filename = tmp_hashes_filename.str();
	      tmp_counts_filename.str("");
	      tmp_counts_filename << seq_dir << "/" << s_org << "/" << "counts.bin";
	      params.tmp_counts_filename = tmp_counts_filename.str();
	      params.num_threads = parallel_for_threads;
	      KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
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

    tp.wait_for_all_tasks();
    for (vector<KMerge::BuilderTask*>::iterator iter=task_ptrs.begin(); iter!=task_ptrs.end(); iter++) {
      delete *iter;
    }
    (void)closedir(dirp);
  }

  delete kmerge;
  

  return 0;
}
