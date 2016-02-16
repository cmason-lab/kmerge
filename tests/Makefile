HOME_DIR = /home/darryl
DEVELOPMENT_INCLUDE_DIR = $(HOME_DIR)/Development/include
DEVELOPMENT_LIB_DIR = $(HOME_DIR)/Development/lib
TEST_INCLUDE_DIR = $(HOME_DIR)/Development/kmerge/tests/src
PROGRAM_SRC_DIR = $(HOME_DIR)/Development/kmerge/src

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -I $(HOME_DIR)/include -I $(DEVELOPMENT_INCLUDE_DIR) -I $(PROGRAM_SRC_DIR)

#OPT = -O2

# Flags passed to the C++ compiler.
#CXXFLAGS += -g -DHAVE_CONFIG_H -Wall -Wextra $(OPT) -DSEQAN_HAS_ZLIB=1 -DDLIB_NO_GUI_SUPPORT -pthread
CXXFLAGS += -g -Wall -Wextra $(OPT) -DDLIB_NO_GUI_SUPPORT -pthread -w

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = kmerge_unittest

#CXX=$(HOME_DIR)/hdf5-1.8.13/bin/h5c++ -std=c++11 -mssse3
CXX=g++ -std=c++11

# libraries to include

LDFLAGS+=-lm -lrt -lm -lz -lcityhash -L $(HOME_DIR)/lib -L $(HOME_DIR)/lib64 -L $(DEVELOPMENT_LIB_DIR) -lpthread

# House-keeping build targets.

all : $(TESTS)

clean :
	rm -f $(TESTS) *.o

kmerge_unittest : kmerge_unittest.o kmerge.o murmur_hash.o spookyV2.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(HOME_DIR)/include/dlib/all/source.cpp $^ -o $@ $(LDFLAGS)

kmerge.o: $(PROGRAM_SRC_DIR)/kmerge.h $(PROGRAM_SRC_DIR)/kmerge.cc $(DEVELOPMENT_INCLUDE_DIR)/lookup3.c $(DEVELOPMENT_INCLUDE_DIR)/hash.c
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(PROGRAM_SRC_DIR)/kmerge.cc -o $@ $(LDFLAGS)

spookyV2.o: $(DEVELOPMENT_INCLUDE_DIR)/SpookyV2.h $(DEVELOPMENT_INCLUDE_DIR)/SpookyV2.cpp	
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(DEVELOPMENT_INCLUDE_DIR)/SpookyV2.cpp -o $@ $(LDFLAGS)

murmur_hash.o: $(DEVELOPMENT_INCLUDE_DIR)/MurmurHash3.h $(DEVELOPMENT_INCLUDE_DIR)/MurmurHash3.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(DEVELOPMENT_INCLUDE_DIR)/MurmurHash3.cpp -o $@ $(LDFLAGS)

kmerge_unittest.o : $(TEST_INCLUDE_DIR)/kmerge_unittest.cc $(DEVELOPMENT_INCLUDE_DIR)/lookup3.c $(PROGRAM_SRC_DIR)/kmerge.h $(PROGRAM_SRC_DIR)/kmerge.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(TEST_INCLUDE_DIR)/kmerge_unittest.cc -o $@ $(LDFLAGS)

