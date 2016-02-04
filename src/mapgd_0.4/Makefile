#Testing with different standards. 
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y
#Right now I need c++11 to support constructor delegation, which I like a lot.

MODE=-O3
STD=-std=c++11
DESTDIR=/usr/local/bin
LOCAL=../../bin/mapgd
WARNINGS=-Wall -Werror=return-type -Wno-write-strings
VERSION=-D VERSION=\"0.4.1\"
CXXFLAGS=$(VERSION) $(STD) $(MODE) $(WARNINGS) -fopenmp -fno-exceptions -pedantic
LDLIBS=-l sqlite3
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
	genotype.o \
	estimate-pooled.o \
	sam2idx.o \
	proview.o \
	interface.o \
	map2genotype.o \
	relatedness.o \
	sql/writesql.o 
	

OBJS=allele_stat.o \
	base.o \
	binomial.o \
	data_types/data.o \
	file-index.o \
	filter.o \
	individual-likelihood.o \
	lnmultinomial.o \
	locus.o \
	newton-method-theta.o \
	map-file.o \
	models.o \
	pooled_data.o \
	pooled-likelihood.o \
	quartet.o \
	relatedness_data.o \
	sample_gof.o \
	data_types/sample_name.o \
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

noomp: CXXFLAGS:=$(NOOMP) $(VERSION) $(STD) $(MODE) $(WARNINGS)
noomp: all

debug: CXXFLAGS:=$(DEBUG) $(VERSION) $(STD) $(MODE) $(WARNINGS)
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

