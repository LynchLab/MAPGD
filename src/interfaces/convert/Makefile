CXXFLAGS=-std=c++11 -Wall -W -O3  

all: dgrid.o genotype.o readtable.o convert.o
	$(CXX) $^ $(CXXFLAGS) -o ../../../bin/convert
clean:
	rm -f *.o ../../../bin/convert
