#include "kmerge.h"
#include "lookup3.c"

using namespace std;

uint hashKmer(const string kMer) {
  return((uint) hashlittle(kMer.c_str(), kMer.length(), 0));
}
