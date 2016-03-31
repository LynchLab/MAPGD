#Testing with different standards. 
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y
#Right now I need c++11 to support constructor delegation, which I like a lot.

MODE=-O3
STD=-std=c++11
DESTDIR=/usr/local/bin
DOCPATH=\"$(shell pwd)/../../docs/man/\"
WARNINGS=-Wall -Werror=return-type -Wno-write-strings
MINOR=`cat VERSION | cut -d '-' -f 1`
VERSION=-D VERSION=\"0.4.$(MINOR)\"
INCLUDE=-I ./ -I ./data_types/ -I ./commands/ -I ./raw/ 

#OPENMP
	
CXXFLAGS=$(VERSION) $(STD) $(MODE) $(WARNINGS) -fno-exceptions -pedantic $(INCLUDE) 

DEBUG=-g -D DEBUG
NOOMP=-D NOOMP
NOSQL=-D NOSQL
NOGSL=-D NOGSL

OBJDIR=obj

ifeq ($(CXX), i686-w64-mingw32-c++)
	TARGET=../../bin/mapgd-win32.exe
	CXXFLAGS+= -static
	HAVE_SQL=false	
	HAVE_OMP=false
	HAVE_GSL=false
else ifeq ($(CXX), x86_64-w64-mingw32-g++)
	TARGET=../../bin/mapgd-win64.exe
	CXXFLAGS+= -static
	HAVE_SQL=false
	HAVE_OMP=false
	HAVE_GSL=false
else ifeq ($(CXX), clang)
	TARGET=../../bin/mapgd
	HAVE_OMP=false
else
	TARGET=../../bin/mapgd
endif

ifeq ($(HAVE_SQL), yes)
	INCLUDE += -I ./sql
	LDLIBS += -l sqlite3 
else
	CXXFLAGS += $(NOSQL)
endif

ifeq ($(HAVE_OMP), yes)
	CXXFLAGS += $(OPENMP_CXXFLAGS)
else
	CXXFLAGS += $(NOOMP)
endif

ifeq ($(HAVE_GSL), yes)
	LDLIBS +=-lgsl -lgslcblas -lm 
else
	CXXFLAGS += $(NOGSL)
endif

ifeq ($(HAVE_MAN), yes)
	CXXFLAGS += -D USE_MAN -D PATH=$(DOCPATH)
endif

DOCDIR=../../docs

SUBDIR=raw sql data_types commands data_conversion
#>! Just a list of flags to try at some later date.
# -oFast -ffloat-store -fexcess-precision=style -ffast-math -fno-rounding-math -fno-signaling-nans
# -fcx-limited-range -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math
# -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math
# -frounding-math -fsingle-precision-constant -fcx-limited-range -fcx-fortran-rules

#>! TODO ORGANIZE THESE!!

MAPGD=mapgd.o

#These all depend on interface.h
COMMAND_OBJ=estimate-individual.o \
	estimate-pooled.o \
	mapgd-help.o \
	sam2idx.o \
	proview.o \
	interface.o \
	map2genotype.o \
	PopLD.o \
	filter.o \
	relatedness.o

#Depend on db_api
SQL_OBJ=sql/writesql.o \
	sql/readsql.o \
	sql/db_api.o 

BASE_OBJ=binomial.o \
	individual-likelihood.o \
	lnmultinomial.o \
	newton-method-theta.o \
	map-file.o \
	models.o \
	pooled-likelihood.o \
	stream-tools.o 

#These all depend on data.h
DATA_OBJ=data_types/allele.o \
	data_types/data.o \
	data_types/file-index.o \
	data_types/linkage_data.o \
	data_types/pooled_data.o \
	data_types/population.o \
	data_types/sample_gof.o \
	data_types/relatedness_data.o \
	data_types/locus.o \
	data_types/sample_name.o 

RAW_OBJ=raw/genotype.o \
	raw/quartet.o \
	raw/genotype_pair.o \
	raw/base.o

CONVERSION_OBJ=data_conversion/vcf_file.o


ifeq ($(HAVE_SQL), yes)
	OBJ=$(MAPGD) $(COMMAND_OBJ) $(SQL_OBJ) $(BASE_OBJ) $(RAW_OBJ) $(DATA_OBJ) $(CONVERSION_OBJ)
else
	OBJ=$(MAPGD) $(COMMAND_OBJ) $(BASE_OBJ) $(RAW_OBJ) $(DATA_OBJ)
endif

all: $(OBJ)
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDLIBS)
	@echo "\n\nConsider typing 'make install DESTDIR=~/bin/' to install a local copy\n"

DATA_SOURCE=data_types/*.cc

data.h: data_types/data.cc
	@echo "updating data_types"
	touch $(DATA_SOURCE)

interface.h: interface.cc
	@echo "updating commands"
	touch $(COMMAND_SOURCE)

noomp: CXXFLAGS:=$(CXXFLAGS) $(NOOMP)
noomp: all

debug: CXXFLAGS:=$(CXXFLAGS) $(DEBUG)
debug: all

increment:
	echo -n $$(($(MINOR)+1)) > VERSION
	echo -n '-' >> VERSION
	date +%F >> VERSION
.PHONY: increment

install: all
	install -m 0755 $(TARGET) $(DESTDIR)
.PHONY: install


docs:
	$(eval COMMAND=`$(TARGET) -a`) 
	echo -n $(COMMAND)
	for COM in $(COMMAND); do \
		if [ ! -f $(DOCDIR)/man/man1/mapgd-$$COM.1.gz ]; then \
			echo $(TARGET) $$COM \$$1 > mapgd-$$COM;  \
			chmod 777 mapgd-$$COM;  \
			help2man ./mapgd-$$COM > $(DOCDIR)/man/man1/mapgd-$$COM.1; \
			gzip $(DOCDIR)/man/man1/mapgd-$$COM.1; \
			rm mapgd-$$COM; \
		fi \
	done
#.PHONY: clean
#	help2man

clean:
	rm -f *.o
	for DIR in $(SUBDIR); do \
		cd $$DIR; \
		rm -f *.o; \
		cd ..; \
	done
.PHONY: clean
