CXXFLAGS= -g -std=c++0x -O3 -Wall -Werror=return-type -Wno-write-strings -fopenmp
prefix=/usr/local

mapdg: mapgd.o allele_stat.o binomial.o compare-pooled.o estimate-individual.o estimate-pooled.o genotype.o individual-likelihood.o interface.o likelihood.o lnmultinomial.o locus.o newton-method-theta.o pooled-likelihood.o population.o pro-file.o proview.o quartet.o stream-tools.o models.o
	$(CXX) $^ $(CXXFLAGS) -o ../../bin/mapgd

estimate: allele_stat.o binomial.o estimate-individual.o estimate-pooled.o genotype.o individual-likelihood.o interface.o likelihood.o lnmultinomial.o locus.o pooled-likelihood.o population.o pro-file.o quartet.o stream-tools.o models.o
	$(CXX) $^ $(CXXFLAGS) -o ../../bin/mapgd-estimate

compare: allele_stat.o binomial.o compare-pooled.o o genotype.o individual-likelihood.o interface.o likelihood.o lnmultinomial.o locus.o newton-method-theta.o pooled-likelihood.o population.o pro-file.o quartet.o stream-tools.o models.o
	$(CXX) $^ $(CXXFLAGS) -o ../../bin/mapgd-compare

proview: mapgd.o allele_stat.o binomial.o compare-pooled.o estimate-individual.o estimate-pooled.o genotype.o individual-likelihood.o interface.o likelihood.o lnmultinomial.o locus.o newton-method-theta.o pooled-likelihood.o population.o pro-file.o proview.o quartet.o stream-tools.o models.o
	$(CXX) $^ $(CXXFLAGS) -o ../../bin/mapgd-proview



.PHONY: debug

.PHONY: debug
debug: CXXFLAGS=-g -std=c++11 -Wall -W -Og
debug: mapdg


install: mapdg
	install -m 0755 ../../bin/mapgd $(prefix)/bin
.PHONY: install
	
clean:
	rm -f *.o
	rm -f ../../bin/*

