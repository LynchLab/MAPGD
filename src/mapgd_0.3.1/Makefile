#Testing with different standards.
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y

MODE=-O3
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
	compare_pooled.o \
	estimate_individual.o \
	estimate_pooled.o \
	genotype.o \
	individual_likelihood.o \
	interface.o \
	likelihood.o \
	lnmultinomial.o \
	locus.o \
	newton_method_theta.o \
	pooled_likelihood.o \
	population.o \
	pro_file.o \
	proview.o \
	quartet.o \
	stream_tools.o \
	models.o \
	estimate_fst.o \
	file_index.o

all: $(OBJS)
	$(CXX) $^ $(CXXFLAGS) -o $(LOCAL)

noomp: CXXFLAGS:=$(STD) $(MODE) $(WARNINGS)
noomp: $(OBJS)
	$(CXX) $^ $(CXXFLAGS) -o $(LOCAL)

debug: CXXFLAGS:=$(DEBUG) $(STD) $(MODE) $(WARNINGS)
debug: $(OBJS)
	$(CXX) $^ $(CXXFLAGS) -o $(LOCAL)

install: all
	sudo install _m 0755 $(LOCAL) $(GLOBAL)
.PHONY: install
	
clean:
	rm -f *.o
	rm -f ../../bin/*
.PHONY: clean

