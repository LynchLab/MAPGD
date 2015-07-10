#CXXFLAGS= -g -std=c++0x -Wall -W -O3 -Wno-write-strings
CXXFLAGS= -g -std=c++0x -O3 -Wno-write-strings -fopenmp
prefix=/usr/local

mapdg: mapgd.o interface.o proview.o proFile.o estimatePooled.o comparePooled.o streamtools.o pooledLikelihood.o estimateInd.o indLikelihood.o Likelihood.o calcInd.o  
	$(CXX) $^ $(CXXFLAGS) -o ../../bin/mapgd

.PHONY: test
test:
	../scripts/RunMeLinux.sh ../test/test.pileup

.PHONY: debug
debug: CXXFLAGS=-g -std=c++11 -Wall -W -Og
debug: mapdg


install: mapdg
	install -m 0755 ../../bin/mapgd $(prefix)/bin
.PHONY: install
	
clean:
	rm -f *.o
	rm -f ../../bin/*

