SRC=src
HOME_DIR=/home/darryl
CXX=/home/darryl/hdf5/bin/h5c++ 
CXX+=-W -Wall -Wno-long-long -pedantic -Wno-variadic-macros
#LDFLAGS+=-L/home/darryl/lib -lz -lrt 
LDFLAGS+=-L/home/darryl/lib -lfastquery_nompi -lm -L/home/darryl/hdf5/lib -lhdf5 -lhdf5_hl -lz -L/home/darryl/fastbit-ibis1.3.8/lib -lfastbit -lrt -L/home/darryl/lib -lsz -L $(HOME_DIR)/fastquery-0.8.2.8/lib
CXXFLAGS += -I $(HOME_DIR)/include -DFQ_NOMPI -I $(HOME_DIR)/hdf5/include -I $(HOME_DIR)/fastbit-ibis1.3.8/win -I $(HOME_DIR)/fastbit-ibis1.3.8/include -I $(HOME_DIR)/fastquery-0.8.2.8/include -D_NOMPI -DFQ_NOMPI

default: all
all: build_kmer_matrix

build_kmer_matrix: main.o
	$(CXX) $(LDFLAGS) -o build_kmer_matrix main.o

main.o: $(SRC)/kmerge.cc $(SRC)/lookup3.c
	$(CXX) $(CXXFLAGS) -c -o main.o $(SRC)/kmerge.cc

clean:
	rm -f main.o main

.PHONY: default all clean