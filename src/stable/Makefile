CXXFLAGS= -g -std=c++0x -O3 -Wno-write-strings -fopenmp
prefix=/usr/local

mapdg: mapgd.o interface.o proview.o proFile.o estimatePooled.o comparePooled.o streamtools.o pooledLikelihood.o estimateInd.o indLikelihood.o Likelihood.o   
	$(CXX) $^ $(CXXFLAGS) -o ../../bin/mapgd

.PHONY: debug
debug: CXXFLAGS=-g -std=c++11 -Wall -W -Og
debug: mapdg


install: mapdg
	install -m 0755 ../../bin/mapgd $(prefix)/bin
.PHONY: install
	
clean:
	rm -f *.o
	rm -f ../../bin/*

