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

using namespace seqan;
using namespace std;


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

  //Parse command line.
  ArgumentParser::ParseResult res = parse(parser, argc, argv);

  // If parsing was not successful then exit with code 1 if there were errors.
  // Otherwise, exit with code 0 (e.g. help was printed).

  if (res != ArgumentParser::PARSE_OK)
    return res == ArgumentParser::PARSE_ERROR;

  uint k_val_start = 0;
  uint k_val_end = 0;
  H5std_string hdf5_file_name;
  getArgumentValue(hdf5_file_name, parser, 0);
  getArgumentValue(k_val_start, parser, 1);
  getArgumentValue(k_val_end, parser, 2);

  
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
  

  dirp = opendir(".");
  KMerge *kmerge = new KMerge(hdf5_file_name);
  while ((dp = readdir(dirp)) != NULL) {
    if(stat(dp->d_name, &st) == 0) {
      if (S_ISDIR(st.st_mode)) {
	string s_org(dp->d_name);
	if (s_org.compare(".") != 0 && s_org.compare("..") != 0) {
	  cout << "Working on " << s_org << endl;
	  hashes.clear();
	  counts.clear();
	  try {
	    for (uint k = k_val_start; k <= k_val_end; k=k+2) {
	      file_name.str("");
	      file_name << "k" << k << ".counts.gz";
	      file_loc.str("");
	      file_loc << s_org << "/" << file_name.str();
	      if (!(kmerge->parseKmerCountsFile(file_loc.str(), hashes, counts))) {
		cerr << "Unable to parse: " << file_loc.str() << endl;
	      } else {
		cout << "Finished parsing: " << file_loc.str() << endl;
		cout << "Hashes vector size now: " << hashes.size() << endl;
	      }
	    }
	    group_name = s_org;
	    dataset_name.str("");
	    dataset_name << "/" << s_org << "/" << "kmer_hash";
	    hash_dataset_name = dataset_name.str();
	    if (!(kmerge->addDatasetToHDF5File(group_name, hash_dataset_name, hashes.size(), &hashes[0], true))) {
	      cerr << "Unable to add hashes for " << s_org << endl;
	    }
	    dataset_name.str("");
	    dataset_name << "/" << s_org << "/" << "count";
	    count_dataset_name = dataset_name.str();
	    if (!(kmerge->addDatasetToHDF5File(group_name, count_dataset_name, counts.size(), &counts[0], false))) {
	      cerr << "Unable to add counts for " << s_org << endl;
	    }
	    nz_count += hashes.size();
	    // done processing kmers for this organism                                             
	    org_count++;
	    cout << hashes.size() << " k-mer hashes for " << s_org << endl;
	    cout << "Done." << endl;
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

  cout << "Org count: " << org_count << endl;
  cout << "Matrix Non-Zeros: " << nz_count << endl;
  delete kmerge;
  (void)closedir(dirp);

  return 0;
}
