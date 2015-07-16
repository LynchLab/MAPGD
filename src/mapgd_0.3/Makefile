CXXFLAGS= -g -std=c++0x -O3 -Wall -Wno-write-strings -fopenmp
prefix=/usr/local

mapdg: mapgd.o stream-tools.o interface.o proview.o estimate-pooled.o compare-pooled.o estimate-individual.o gcf-file.o gcf-file-header.o pro-file.o quartet.o locus.o likelihood.o pooled-likelihood.o individual-likelihood.o genotype.o population.o lnmultinomial.o binomial.o allele_stat.o newton-method-theta.o hypergeometric_pFq.o
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

