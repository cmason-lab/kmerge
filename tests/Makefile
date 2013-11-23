HOME_DIR = /home/darryl
DEVELOPMENT_INCLUDE_DIR = $(HOME_DIR)/Development/include
TEST_INCLUDE_DIR = $(HOME_DIR)/Development/kmerge/tests/src
PROGRAM_SRC_DIR = $(HOME_DIR)/Development/kmerge/src

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -I $(HOME_DIR)/fastbit-ibis1.3.8/win -I $(HOME_DIR)/fastbit-ibis1.3.8/include -D_NOMPI -DFQ_NOMPI -I $(HOME_DIR)/include -I -DFQ_NOMPI -I $(HOME_DIR)/fastquery-0.8.2.8/include -I $(DEVELOPMENT_INCLUDE_DIR) -I $(PROGRAM_SRC_DIR)

#OPT = -O2

# Flags passed to the C++ compiler.
CXXFLAGS += -g -DHAVE_CONFIG_H -Wall -Wextra $(OPT) 

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = kmerge_unittest

# using hdf5 tools for compiling
CXX=/home/darryl/hdf5/bin/h5c++
#CXX=g++

# libraries to include

LDFLAGS+=-L $(HOME_DIR)/fastquery-0.8.2.8/lib -lfastquery_nompi -lm -L $(HOME_DIR)/fastbit-ibis1.3.8/lib -lfastbit -lrt -lm -lz -L $(HOME_DIR)/lib -lsz -L $(HOME_DIR)/lib64 -larmadillo

# House-keeping build targets.

all : $(TESTS)

clean :
	rm -f $(TESTS) *.o

kmerge_unittest : kmerge_unittest.o kmerge.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

kmerge.o: $(PROGRAM_SRC_DIR)/kmerge.h $(DEVELOPMENT_INCLUDE_DIR)/lookup3.c
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(PROGRAM_SRC_DIR)/kmerge.cc -o $@

kmerge_unittest.o : $(TEST_INCLUDE_DIR)/kmerge_unittest.cc $(DEVELOPMENT_INCLUDE_DIR)/lookup3.c $(PROGRAM_SRC_DIR)/kmerge.h $(PROGRAM_SRC_DIR)/kmerge.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(TEST_INCLUDE_DIR)/kmerge_unittest.cc $(LDFLAGS)
