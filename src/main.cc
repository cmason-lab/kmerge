#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/arg_parse.h>
#include <dirent.h>
#include <sys/stat.h>
#include <string>
#include <seqan/stream.h>
#include <fstream>
#include "kmerge.h"
#include <dlib/threads.h>

using namespace seqan;
using namespace std;
using namespace dlib;

int main(int argc, char const ** argv) {
  ArgumentParser parser("kmerge");


  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "STRING"));
  setHelpText(parser, 0, "HDF5 file name");
  addArgument(parser, ArgParseArgument(ArgParseArgument::INTEGER, "INT"));
  setHelpText(parser, 1, "Start k-mer value >= 5");
  addArgument(parser, ArgParseArgument(ArgParseArgument::INTEGER, "INT"));
  setHelpText(parser, 2, "End k-mer value <= 31");
  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "STRING"));
  setHelpText(parser, 3, "Full path to location containing directories of sequences and taxonomy information");
  addOption(parser, ArgParseOption("t", "threads", "Number of threads to use.",
					  ArgParseArgument::INTEGER, "INT"));
  addOption(parser, ArgParseOption("h", "hash_func", "Hash function to use for k-mers",
				   ArgParseArgument::INTEGER, "STRING"));
  setMinValue(parser, "t", "1");
  setMaxValue(parser, "t", "80");

  //Parse command line.
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  // If parsing was not successful then exit with code 1 if there were errors.
  // Otherwise, exit with code 0 (e.g. help was printed).

  if (res != ArgumentParser::PARSE_OK)
    return res == ArgumentParser::PARSE_ERROR;

  uint k_val_start = 0;
  uint k_val_end = 0;
  std::string hdf5_filename, seq_dir, hash_func = "lookup3";
  uint num_threads = 1;

  getArgumentValue(hdf5_filename, parser, 0);
  getArgumentValue(k_val_start, parser, 1);
  getArgumentValue(k_val_end, parser, 2);
  getArgumentValue(seq_dir, parser, 3);
  if (isSet(parser, "t")) {
    getOptionValue(num_threads, parser, "t");
  }
  if (isSet(parser, "h")) {
    getOptionValue(hash_func, parser, "h");
  }

  if ((k_val_start < 5) || (k_val_end > 31)) {
    std::cerr << "k-mer range is out of bounds" << std::endl << "k must be between 5 and 31 (inclusive)" << std::endl;
    return 1;
  }
  
  struct stat st;
  DIR *dirp;
  struct dirent *dp;
  stringstream seq_filename, file_loc, dataset_name;
  vector<uint> hashes;
  vector<uint> counts;
  std::string group_name("");
  std::string hash_dataset_name("");
  std::string count_dataset_name("");
  vector<KMerge::BuilderTask*> task_ptrs;
 
  thread_pool tp(num_threads);
  dirp = opendir(seq_dir.c_str());
  KMerge *kmerge = new KMerge(hdf5_filename, hash_func);
  while ((dp = readdir(dirp)) != NULL) {
    if(stat(dp->d_name, &st) == 0) {
      if (S_ISDIR(st.st_mode)) {
	string s_org(dp->d_name);
	if (s_org.compare(".") != 0 && s_org.compare("..") != 0) {
	  try {
	    param_struct params;
	    params.kmerge = kmerge;
	    params.hdf5_filename = hdf5_filename;
	    params.k_val_start = k_val_start;
	    params.k_val_end = k_val_end;
	    params.group_name = s_org;
	    seq_filename.str("");
	    seq_filename << seq_dir << "/" << s_org << "/" << s_org << ".fasta.gz" << endl;
	    params.seq_filename = seq_filename.str();
	    dataset_name.str("");
	    dataset_name << "/" << s_org << "/" << "kmer_hash";
	    params.hash_dataset_name = dataset_name.str();
	    dataset_name.str("");
	    dataset_name << "/" << s_org << "/" << "count";
	    params.counts_dataset_name = dataset_name.str();
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

  delete kmerge;
  
  (void)closedir(dirp);

  return 0;
}
