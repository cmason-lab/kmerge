#include "H5Cpp.h"
#include "queryProcessor.h"
#include "gtest/gtest.h"

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
using std::cout;
using std::endl;
using std::string;
#endif  // H5_NO_STD                                                                                                                                                          
#endif

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif


int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
  //std::cout << "Hello World";
}
