SRC=/home/darryl/Development/kmerge/src
HOME_DIR=/home/darryl
DEVELOPMENT_INCLUDE_DIR=$(HOME_DIR)/Development/include
DEVELOPMENT_LIB_DIR = $(HOME_DIR)/Development/lib
CXX=g++
CXXFLAGS+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -pthread -DDLIB_NO_GUI_SUPPORT -fno-strict-aliasing -DLINUX_ -DOC_NEW_STYLE_INCLUDES -Wno-deprecated -D_REENTRANT

CPPFLAGS += -I $(HOME_DIR)/include -I $(DEVELOPMENT_INCLUDE_DIR) -I $(DEVELOPMENT_INCLUDE_DIR)/PicklingTools151Release/C++/opencontainers_1_7_7/include -I $(DEVELOPMENT_INCLUDE_DIR)/PicklingTools151Release/C++

LDFLAGS+=-lrt -lm -lz -L $(HOME_DIR)/lib -L $(DEVELOPMENT_LIB_DIR) -lptools


default: all
all: pickle

pickle: main.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) main.o $(HOME_DIR)/include/dlib/all/source.cpp -o $@ $(LDFLAGS)

main.o: main.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c main.cc -o $@ $(LDFLAGS)


clean:
	rm -f *.o combine

.PHONY: default all clean