#include <string>
using namespace std;
typedef unsigned int uint;

uint hashKmer(const std::string);
bool createDataset(void);
bool appendToDataset(void);
bool updateDataset(void);
uint * getDatasetValues(void);
