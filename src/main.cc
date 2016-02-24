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

  parser.add_option("k", "Start and end k-mer values", 2);
  parser.add_option("d", "Location of sequences and taxonomy directories", 1);
  parser.add_option("i", "Location of single sequence file", 1);
  parser.add_option("p", "Sequencing reads are paired-end", 0);
  parser.add_option("o", "Location to write output files (must already exist)", 1);
  parser.add_option("u", "Truncated hash value in bits (8 or 16)", 1);
  parser.add_option("t", "Max number of threads to use for worker and hash threads", 2);
  parser.add_option("f", "Hash function to use for k-mers", 1);
  parser.add_option("h","Display this help message.");


  //Parse command line.
  parser.parse(argc,argv);

  if (!(parser.parsed_line())) {
    std::cerr << "Line not parsed" << std::endl;
    return 0;
  }
  //Ensure options only defined once
  const char* one_time_opts[] = {"k", "d", "i", "p", "o", "u", "t", "f", "h"};
  parser.check_one_time_options(one_time_opts);

  //Check only one of d and i options given
  parser.check_incompatible_options("d", "i");
  parser.check_incompatible_options("d", "p");

  parser.check_option_arg_range("k", 3, 31);

  // check if the -h option was given on the command line
  if (parser.option("h")) {
    // display all the command line options
    std::cout << "Usage: kmerge -k k_start k_end (-d|-i|-p|-o|-u|-t|-f)" << std::endl;
    parser.print_options(); 
    return 0;
  }

  // make sure two values input to k
  if (parser.option("k").number_of_arguments() != 2) {
    std::cerr << "Error in command line:\n   You must specify start and end k-mer values for parsing" << std::endl;
    std::cerr << "\nTry the -h option for more information." << std::endl;
    return 0;
  }
 

  uint k_val_start = 0;
  uint k_val_end = 0;
  std::string seq_dir("./"), hash_func("lookup3"), in_file, out_dir("./reference");
  uint work_threads = 1, hash_threads = 1, trunc_mode=0;
  bool paired_end=false;

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
    out_dir = "./sample";
    if (parser.option("p")) {
      paired_end = true;
    }
  } 
  if (parser.option("o")) {
    out_dir = parser.option("o").argument();
  }
  if (parser.option("u")) {
    trunc_mode = atoi(parser.option("u").argument().c_str());
    if (trunc_mode != 8 and trunc_mode != 16) {
      std::cerr << "Error in command line:\n   Hash truncation value must be 8 or 16." << std::endl;
    }
  }
  if (parser.option("t")) {
    // make sure two values input to t
    if (parser.option("t").number_of_arguments() != 2) {
      std::cerr << "Error in command line:\n   You must specify worker and hash thread values if specifying threads" << std::endl;
      std::cerr << "\nTry the -h option for more information." << std::endl;
      return 0;
    }
    work_threads = atoi(parser.option("t").argument(0).c_str());
    if (work_threads < 1) {
      work_threads = 1;
    }
    hash_threads = atoi(parser.option("t").argument(1).c_str());
    if (hash_threads < 1) {
      hash_threads = 1;
    }
  }
  if (parser.option("f")) {
    std::string option = parser.option("f").argument();
    // check that values are valid or use default
    if ((option == "lookup3") || (option == "spooky") || (option == "murmur") || (option == "city")) {
      hash_func = option;
    }
  }

  KMerge *kmerge;

  std::vector<KMerge::BuilderTask*> task_ptrs;
  dlib::thread_pool tp(work_threads);

  if (parser.option("i")) {
    try {
      kmerge = new KMerge(hash_func, "", out_dir, trunc_mode);
      param_struct params;
      params.kmerge = kmerge;
      params.k_val_start = k_val_start;
      params.k_val_end = k_val_end;
      params.group_name = "sample";
      params.seq_filename = in_file;
      params.hashes_filename = out_dir + std::string("/") + params.group_name + std::string(".hashes.bin");
      params.counts_filename = out_dir + std::string("/") + params.group_name + std::string(".counts.bin");
      params.indices_filename = out_dir + std::string("/") + params.group_name + std::string(".indices.bin");
      params.ids_filename = out_dir + std::string("/") + params.group_name + std::string(".ids.bin");
      params.num_threads = hash_threads;
      params.is_ref = false;
      params.paired_end = paired_end;
      KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
      tp.add_task(*task, &KMerge::BuilderTask::execute);
      task_ptrs.push_back(task);
    } catch( exception &e ) {
      cout << e.what() << std::endl;
      return 1;
    }

  } else {
    struct stat st;
    DIR *dirp;
    struct dirent *dp;
    dirp = opendir(seq_dir.c_str());
    kmerge = new KMerge(hash_func, seq_dir, out_dir, trunc_mode);
    while ((dp = readdir(dirp)) != NULL) {
      if(stat(dp->d_name, &st) == 0) {
	if (S_ISDIR(st.st_mode)) {
	  std::string s_org(dp->d_name);
	  if (s_org.compare(".") != 0 && s_org.compare("..") != 0) {
	    try {
	      param_struct params;
	      params.kmerge = kmerge;
	      params.k_val_start = k_val_start;
	      params.k_val_end = k_val_end;
	      params.group_name = s_org;
	      params.seq_filename = seq_dir + params.group_name + std::string("/") + params.group_name + std::string(".fasta.gz");
	      params.hashes_filename = out_dir + std::string("/") + params.group_name + std::string(".hashes.bin");
	      params.counts_filename = out_dir + std::string("/") + params.group_name + std::string(".counts.bin");
	      params.num_threads = hash_threads;
	      params.is_ref = true;
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

    (void)closedir(dirp);
  }

  tp.wait_for_all_tasks();
  for (vector<KMerge::BuilderTask*>::iterator iter=task_ptrs.begin(); iter!=task_ptrs.end(); iter++) {
    delete *iter;
  }

  delete kmerge;
  

  return 0;
}
