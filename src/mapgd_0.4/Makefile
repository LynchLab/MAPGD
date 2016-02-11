#Testing with different standards. 
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y
#Right now I need c++11 to support constructor delegation, which I like a lot.

#MODE=-O3
STD=-std=c++11
DESTDIR=/usr/local/bin
LOCAL=../../bin/mapgd
WARNINGS=-Wall -Werror=return-type -Wno-write-strings
VERSION=-D VERSION=\"0.4.1\"
POSIX=-D POSIX
INCLUDE=-I ./ -I ./data_types/ -I ./commands/ -I ./raw/ -I ./sql
CXXFLAGS=$(VERSION) $(POSIX) $(STD) $(MODE) $(WARNINGS) -fopenmp -fno-exceptions -pedantic $(INCLUDE)
LDLIBS=-l sqlite3 -lgsl -lgslcblas -lm
DEBUG=-g -D DEBUG
NOOMP=-D NOOMP
OBJDIR=obj
CXX=g++

#>! Just a list of flags to try at some later date.

# -oFast -ffloat-store -fexcess-precision=style -ffast-math -fno-rounding-math -fno-signaling-nans
# -fcx-limited-range -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math
# -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math
# -frounding-math -fsingle-precision-constant -fcx-limited-range -fcx-fortran-rules

#>! TODO ORGANIZE THESE!!

COMMANDS=mapgd.o \
	estimate-individual.o \
	estimate-pooled.o \
	sam2idx.o \
	proview.o \
	interface.o \
	map2genotype.o \
	relatedness.o \
	sql/writesql.o 	

OBJS=	base.o \
	binomial.o \
	data_types/allele.o \
	data_types/data.o \
	data_types/file-index.o \
	data_types/linkage_data.o \
	data_types/pooled_data.o \
	data_types/population.o \
	data_types/sample_gof.o \
	data_types/relatedness_data.o \
	data_types/locus.o \
	data_types/sample_name.o \
	filter.o \
	individual-likelihood.o \
	lnmultinomial.o \
	newton-method-theta.o \
	map-file.o \
	models.o \
	pooled-likelihood.o \
	PopLD.o \
	raw/genotype.o \
	raw/quartet.o \
	raw/genotype_pair.o \
	stream-tools.o \
	sql/db_api.o 

#	proview.o \
	compare-pooled.o \
	readtable.o \
	convert.o \
	PopLD.o \

OBJS+=$(COMMANDS)

all: $(OBJS)  
	$(CXX) $^ $(CXXFLAGS) -o $(LOCAL) $(LDLIBS)
	@echo "\n\nConsider typing 'make install DESTDIR=~/bin/' to install a local copy\n"

noomp: CXXFLAGS:=$(CXXFLAGS) $(NOOMP)
noomp: all

debug: CXXFLAGS:=$(CXXFLAGS) $(DEBUG)
debug: all

install: all
	install -m 0755 $(LOCAL) $(DESTDIR)
#	install ? $(MANDIR)/man1/
.PHONY: install
	
clean:
	rm -f *.o
	rm -f data_types/*.o
	rm -f sql/*.o
	rm -f ../../bin/*
	rm -f ~/bin/mapgd
.PHONY: clean

