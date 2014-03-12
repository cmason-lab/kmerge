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
#include "H5Cpp.h"
#include <dlib/threads.h>

using namespace seqan;
using namespace std;
using namespace dlib;

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

int main(int argc, char const ** argv) {
  ArgumentParser parser("kmerge");


  addArgument(parser, ArgParseArgument(ArgParseArgument::STRING, "STRING"));
  setHelpText(parser, 0, "HDF5 file name");
  addArgument(parser, ArgParseArgument(ArgParseArgument::INTEGER, "INT"));
  setHelpText(parser, 1, "Start k-mer value");
  addArgument(parser, ArgParseArgument(ArgParseArgument::INTEGER, "INT"));
  setHelpText(parser, 2, "End k-mer value");
  addOption(parser, ArgParseOption("t", "threads", "Number of threads to use.",
					  ArgParseArgument::INTEGER, "INT"));
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
  H5std_string hdf5_file_name;
  uint num_threads = 1;
  getArgumentValue(hdf5_file_name, parser, 0);
  getArgumentValue(k_val_start, parser, 1);
  getArgumentValue(k_val_end, parser, 2);
  if (isSet(parser, "t")) {
    getOptionValue(num_threads, parser, "t");
  }
  
  struct stat st;
  DIR *dirp;
  struct dirent *dp;
  stringstream file_name, file_loc, dataset_name;
  vector<uint> hashes;
  vector<uint> counts;
  uint org_count = 0;
  uint nz_count = 0;
  H5std_string group_name("");
  H5std_string hash_dataset_name("");
  H5std_string count_dataset_name("");
  //vector<param_struct*> param_ptrs;
  vector<KMerge::BuilderTask*> task_ptrs;
 
  thread_pool tp(num_threads);
  dirp = opendir(".");
  KMerge *kmerge = new KMerge(hdf5_file_name);
  while ((dp = readdir(dirp)) != NULL) {
    if(stat(dp->d_name, &st) == 0) {
      if (S_ISDIR(st.st_mode)) {
	string s_org(dp->d_name);
	if (s_org.compare(".") != 0 && s_org.compare("..") != 0) {
	  try {
	    param_struct params;
	    params.kmerge = kmerge;
	    params.hdf5_file_name = hdf5_file_name;
	    params.k_val_start = k_val_start;
	    params.k_val_end = k_val_end;
	    params.group_name = s_org;
	    dataset_name.str("");
	    dataset_name << "/" << s_org << "/" << "kmer_hash";
	    params.hash_dataset_name = dataset_name.str();
	    dataset_name.str("");
	    dataset_name << "/" << s_org << "/" << "count";
	    params.count_dataset_name = dataset_name.str();
	    KMerge::BuilderTask* task = new KMerge::BuilderTask(params);
	    //param_ptrs.push_back(params);
            tp.add_task(*task, &KMerge::BuilderTask::execute);
	    task_ptrs.push_back(task);
	    //nz_count += hashes.size();
	    // done processing kmers for this organism                                             
	    //org_count++;
	  }
	  catch( FileIException error ) {
	    error.printError();
	    continue;
	  }
	  // catch failure caused by the DataSet operations  
	  catch( DataSetIException error ) {
	    error.printError();
	    continue;
	  }
	  // catch failure caused by the DataSpace operations 
	  catch( DataSpaceIException error ) {
	    error.printError();
	    continue;
	  }
	}
      }
    }
  }
  tp.wait_for_all_tasks();
  
  //cout << "Org count: " << org_count << endl;
  //cout << "Matrix Non-Zeros: " << nz_count << endl;
  for (vector<KMerge::BuilderTask*>::iterator iter=task_ptrs.begin(); iter!=task_ptrs.end(); iter++) {
    delete *iter;
  }

  /*for (vector<param_struct*>::iterator iter=param_ptrs.begin(); iter!=param_ptrs.end(); iter++) {
    delete *iter;
    }*/

  delete kmerge;
  
  (void)closedir(dirp);

  return 0;
}
