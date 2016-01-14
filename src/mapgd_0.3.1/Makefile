#Testing with different standards.
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y

#MODE=-O3
STD=-std=c++0x
GLOBAL=/usr/local/bin
LOCAL=../../bin/mapgd
WARNINGS=-Wall -Werror=return-type -Wno-write-strings
CXXFLAGS=$(STD) $(MODE) $(WARNINGS) -fopenmp
DEBUG= -g 
CXX=g++


#>! TODO ORGANIZE THESE!!
OBJS=	mapgd.o \
	allele_stat.o \
	binomial.o \
	compare-pooled.o \
	estimate-individual.o \
	estimate-pooled.o \
	genotype.o \
	individual-likelihood.o \
	interface.o \
	likelihood.o \
	lnmultinomial.o \
	locus.o \
	newton-method-theta.o \
	pooled-likelihood.o \
	population.o \
	pro-file.o \
	proview.o \
	quartet.o \
	stream-tools.o \
	models.o \
	estimate-fst.o \
	PopLD.o \
	relatedness.o \
	file-index.o

all: $(OBJS)
	$(CXX) $^ $(CXXFLAGS) -o $(LOCAL)

noomp: CXXFLAGS:=$(STD) $(MODE) $(WARNINGS)
noomp: $(OBJS)
	$(CXX) $^ $(CXXFLAGS) -o $(LOCAL)

debug: CXXFLAGS:=$(DEBUG) $(STD) $(MODE) $(WARNINGS)
debug: $(OBJS)
	$(CXX) $^ $(CXXFLAGS) -o $(LOCAL)

install: all
	sudo install -m 0755 $(LOCAL) $(GLOBAL)
.PHONY: install
	
clean:
	rm -f *.o
	rm -f ../../bin/*
.PHONY: clean

