#Testing with different standards.
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y

MODE=-O3
STD=-std=c++0x
GLOBAL=/usr/local/bin
LOCAL=../../bin/mapgd
WARNINGS=-Wall -Werror=return-type -Wno-write-strings
VERSION=-D VERSION=\"0.4.1\"
CXXFLAGS=$(VERSION) $(STD) $(MODE) $(WARNINGS) -fopenmp
DEBUG=-g -D DEBUG
NOOMP=-D NOOMP
OBJDIR=obj
CXX=g++


#>! TODO ORGANIZE THESE!!

OBJS=allele_stat.o \
	base.o \
	binomial.o \
	clone_gof.o \
	convert.o \
	compare-pooled.o \
	estimate-individual.o \
	estimate-pooled.o \
	file-index.o \
	filter.o \
	genotype.o \
	individual-likelihood.o \
	interface.o \
	lnmultinomial.o \
	locus.o \
	newton-method-theta.o \
	models.o \
	map2genotype.o \
	pooled-likelihood.o \
	proview-new.o \
	pro-file.o \
	PopLD.o \
	quartet.o \
	stream-tools.o \
	readtable.o \
	convert.o \
	map-file.o \
	mapgd.o

#	proview.o \

all: $(OBJS)
	$(CXX) $^ $(CXXFLAGS) -o $(LOCAL) 

noomp: CXXFLAGS:=$(NOOMP) $(STD) $(MODE) $(WARNINGS)
noomp: all

debug: CXXFLAGS:=$(DEBUG) $(VERSION) $(STD) $(MODE) $(WARNINGS)
debug: all

install: all
	sudo install -m 0755 $(LOCAL) $(GLOBAL)
.PHONY: install
	
clean:
	rm -f *.o
	rm -f ../../bin/*
.PHONY: clean

