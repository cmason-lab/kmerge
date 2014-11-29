#include <vector> 
#include <inttypes.h>
#include <dlib/serialize.h>
#include <dlib/compress_stream.h>
#include <numeric>
#include <fstream>

std::vector<uint32_t> load_data(std::string bin_filename, bool sum) { 
  std::vector<uint32_t> vec_in;
  dlib::compress_stream::kernel_3b decompressor;
  std::ifstream ifs;
  std::stringstream ss;

  ifs.open (bin_filename.c_str(), std::ios::binary);
  decompressor.decompress(ifs, ss);
  ifs.close();
  dlib::deserialize(vec_in, ss);
  ss.str(std::string());
  if (sum) std::partial_sum(vec_in.begin(), vec_in.end(), vec_in.begin());

  return vec_in; 
} 
