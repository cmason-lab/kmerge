#include <dlib/serialize.h>
#include <dlib/cmd_line_parser.h>
#include <dlib/compress_stream.h>
#include <dlib/logger.h>
#include <dlib/threads.h>
#include <dirent.h>
#include <sys/stat.h>
#include <math.h>
#include <fstream>
#include <numeric>
#include "chooseser.h"

typedef dlib::compress_stream::kernel_3b cs;
#define MAX_UINT_VAL 4294967295 //2^32-1

using namespace std;

class Task {
public:
  std::string s_org;
  dlib::logger dlog;

  void execute() {

    std::vector<uint> hashes_in, counts_in;
    cs compressor, decompressor;
    std::ifstream ifs;
    std::ofstream ofs;
    std::stringstream ss;

    dlog << dlib::LINFO << "Processing";
    dlog << dlib::LINFO << "Reading hashes";
    ifs.open (std::string(s_org + "/hashes.bin").c_str(), ios::binary);
    dlog << dlib::LINFO << "Decompressing hashes";
    decompressor.decompress(ifs, ss);
    ifs.close();
    dlog << dlib::LINFO << "Deserializing hashes";
    dlib::deserialize(hashes_in, ss);
    ss.str(std::string());
    std::partial_sum(hashes_in.begin(), hashes_in.end(), hashes_in.begin());

    dlog << dlib::LINFO << "Pickling hashes";
    Array<uint> p_hashes(MAX_UINT_VAL);

    for (std::vector<uint>::const_iterator it = hashes_in.begin(); it != hashes_in.end(); it++) {
      p_hashes.append(*it);
    }

    hashes_in.clear();
    std::vector<uint>().swap(hashes_in);

    DumpValToFile(p_hashes, std::string(s_org + "/hashes.pkl"), SERIALIZE_P2);
    
    p_hashes.clear();
    Array<uint>().swap(p_hashes);

    dlog << dlib::LINFO << "Finished pickling hashes";

    dlog << dlib::LINFO << "Reading counts";
    ifs.open (std::string(s_org + "/counts.bin").c_str(), ios::binary);
    dlog << dlib::LINFO << "Decompressing counts";
    decompressor.decompress(ifs, ss);
    ifs.close();
    dlog << dlib::LINFO << "Deserializing counts";
    dlib::deserialize(counts_in, ss);
    ss.str(std::string());

    dlog << dlib::LINFO << "Pickling counts";    
    Array<uint> p_counts(MAX_UINT_VAL);

    for (std::vector<uint>::const_iterator it = counts_in.begin(); it != counts_in.end(); it++) {
      p_counts.append(*it);
    }

    counts_in.clear();
    std::vector<uint>().swap(counts_in);

    DumpValToFile(p_counts, std::string(s_org + "/counts.pkl"), SERIALIZE_P2);
    
    p_counts.clear();
    Array<uint>().swap(p_counts);

    dlog << dlib::LINFO << "Finished pickling counts";
    
    dlog << dlib::LINFO << "Finished processing";
    
  }
    
  Task(const std::string& s_org) : s_org(s_org), dlog(s_org.c_str()) {
    dlog.set_level(dlib::LALL);
  }

  ~Task() {
  }
};


int main(int argc, char const ** argv) {
  dlib::command_line_parser parser;

  parser.add_option("d", "Directory where hashes and count groups are located", 1);
  parser.add_option("t", "Number of threads to use", 1);

  //Parse command line.
  parser.parse(argc,argv);

  if (!(parser.parsed_line())) {
    std::cerr << "Line not parsed" << std::endl;
    return 0;
  }

  std::string dir("./");
  uint num_threads = 1;

  if (parser.option("d")) {
    dir = parser.option("d").argument();
  }
  if (parser.option("t")) {
    num_threads = atoi(parser.option("t").argument().c_str());
  }

  dlib::thread_pool tp(num_threads);
  std::vector<Task*> task_ptrs;
  struct stat st;
  DIR *dirp;
  struct dirent *dp;
  dirp = opendir(dir.c_str());
  while ((dp = readdir(dirp)) != NULL) {
    if(stat(dp->d_name, &st) == 0) {
      if (S_ISDIR(st.st_mode)) {
	std::string s_org(dp->d_name);
	if (s_org.compare(".") != 0 && s_org.compare("..") != 0) {
	  Task* task = new Task(s_org);
	  tp.add_task(*task, &Task::execute);
	  task_ptrs.push_back(task);
	}
      }
    }
  }

  tp.wait_for_all_tasks();
  for (std::vector<Task*>::iterator iter=task_ptrs.begin(); iter!=task_ptrs.end(); iter++) {
    delete *iter;
  }


  return 0;
}
