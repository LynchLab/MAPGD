#Testing with different standards. 
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y
#Right now I need c++11 to support constructor delegation, which I like a lot.

#MODE=-O3
#STD=-std=c++11

WARNINGS=-Wno-write-strings
#WARNINGS=-Wall -Werror=return-type -Wno-write-strings
MINOR=`cat VERSION | cut -d - -f 1`
VERSION=-D VERSION=\"0.4.$(MINOR)\"
INCLUDE=-I ./ -I ./data_types/ -I ./commands/ -I ./raw/ -I ./data_conversion/ -I ./gzstream/ -I ./stream_tools/ -I ./io/ -I ./math_functions/ 

CXXFLAGS += $(VERSION) $(MODE) $(WARNINGS) -fno-exceptions -pedantic $(INCLUDE) -MMD
LDLIBS += -lgzstream -lz -L ./gzstream/ 

DEBUG=-pg -D DEBUG -ggdb
NOOMP=-D NOOMP
NOSQL=-D NOSQL
NOGSL=-D NOGSL
NOHTS=-D NOHTS
NOMPI=-D NOMPI
NOLZ4=-D NOLZ4

#HTS_LIB_PATH=/usr/include/htslib/

OBJDIR=objs

SUBDIR=raw sql data_types commands data_conversion stream_tools io math_functions

VPATH = $(SUBDIR)

ifeq ($(CXX), i686-w64-mingw32-g++)
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
	CXXFLAGS += -D USE_MAN -D PATH=\"$(DOCPATH)\"
endif

ifeq ($(HAVE_HTS), yes)
	INCLUDE += $(HTSLIB_CPPFLAGS)/
	INCLUDE += $(HTSLIB_LDFLAGS)
	LDLIBS += -lhts
else
	CXXFLAGS += $(NOHTS)
endif

ifeq ($(HAVE_LZ4), yes)
	INCLUDE += $(LZ4_INCLUDE)
	#-I ~/src/lz4/lib
	LDLIBS += $(LZ4_LIB)
	#-L ~/src/lz4/lib -llz4
else
	CXXFLAGS += $(NOLZ4)
endif

DOCDIR=../../docs


ifeq ($(UNSTABLE), yes)
	SUBDIR += testing
endif
#>! Just a list of flags to try at some later date.
# -oFast -ffloat-store -fexcess-precision=style -ffast-math -fno-rounding-math -fno-signaling-nans
# -fcx-limited-range -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math
# -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math
# -frounding-math -fsingle-precision-constant -fcx-limited-range -fcx-fortran-rules

#>! TODO ORGANIZE THESE!!

MAPGD=mapgd.o

#These all depend on interface.h

#Headers with templates or inline code

COMMAND_OBJ=$(OBJDIR)/estimate_individual.o \
	$(OBJDIR)/estimate_pooled.o \
	$(OBJDIR)/mapgd_help.o \
	$(OBJDIR)/sam2idx.o \
	$(OBJDIR)/proview.o \
	$(OBJDIR)/interface.o \
	$(OBJDIR)/map2genotype.o \
	$(OBJDIR)/linkage_disequilibrium.o \
	$(OBJDIR)/filter.o \
	$(OBJDIR)/filter_genotype.o \
	$(OBJDIR)/filter_pool.o \
	$(OBJDIR)/fastview.o \
	$(OBJDIR)/relatedness.o \
	$(OBJDIR)/test_keys.o \
	$(OBJDIR)/reml.o \
	$(OBJDIR)/relatedness_test.o


#Depend on htslib 
HTS_OBJ=$(OBJDIR)/vcf-file.o \
	$(OBJDIR)/read_vcf.o \
	$(OBJDIR)/make-vcf2.o \
	$(OBJDIR)/make-vcf.o

#Depend on lz4
LZ4_OBJ=$(OBJDIR)/state.o \
	$(OBJDIR)/simulate.o 

#Depend on db_api
SQL_OBJ=$(OBJDIR)/writesql.o \
	$(OBJDIR)/readsql.o \
	$(OBJDIR)/db_api.o 

#many objects depend on these
BASE_OBJ=$(OBJDIR)/binomial.o \
	$(OBJDIR)/individual_likelihood.o \
	$(OBJDIR)/lnmultinomial.o \
	$(OBJDIR)/newton-method-theta.o \
	$(OBJDIR)/newton-method-ld.o \
	$(OBJDIR)/map_file.o \
	$(OBJDIR)/models.o \
	$(OBJDIR)/pooled_likelihood.o \
	$(OBJDIR)/tmp_buffer.o \
	$(OBJDIR)/stream_tools.o 

#These all depend on data.h
DATA_OBJ=$(OBJDIR)/data.o \
	$(OBJDIR)/allele.o \
	$(OBJDIR)/file_index.o \
	$(OBJDIR)/linkage_data.o \
	$(OBJDIR)/pooled_data.o \
	$(OBJDIR)/population.o \
	$(OBJDIR)/sample_gof.o \
	$(OBJDIR)/relatedness_data.o \
	$(OBJDIR)/bcf2pro-file.o \
	$(OBJDIR)/locus.o \
	$(OBJDIR)/key.o \
	$(OBJDIR)/sample_name.o 

#??
RAW_OBJ=$(OBJDIR)/genotype.o \
	$(OBJDIR)/quartet.o \
	$(OBJDIR)/genotype_pair.o \
	$(OBJDIR)/region.o \
	$(OBJDIR)/constants.o \
	$(OBJDIR)/base.o


OBJ=$(MAPGD) $(COMMAND_OBJ) $(BASE_OBJ) $(RAW_OBJ) $(DATA_OBJ) $(CONVERSION_OBJ) 

ifeq ($(HAVE_MPI), yes)
	MPI_OBJ=$(OBJDIR)/mpi_relatedness.o
	OBJ := $(OBJ) $(MPI_OBJ)
	MPI_LDLIBS += -lmpi -pthread
	FINAL_CXX=$(MPICXX)
else
	MPI_OBJ =
	MPI_LDLIBS =
	FINAL_CXX = $(CXX)
	CXXFLAGS += $(NOMPI)
endif

ifeq ($(HAVE_SQL), yes)
	OBJ := $(OBJ) $(SQL_OBJ)
endif

ifeq ($(HAVE_HTS), yes)
	OBJ := $(OBJ) $(HTS_OBJ)
endif

ifeq ($(HAVE_LZ4), yes)
	OBJ := $(OBJ) $(LZ4_OBJ)
endif

all : gzstream $(TARGET)
#	cd python_api
#	python setup.py clean
#	python setup.py install

$(OBJDIR)/%.o: %.cc
	$(FINAL_CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

#$(TARGET) 

gzstream :
	$(MAKE) -C gzstream 
.PHONY: gzstream

$(TARGET) : $(OBJ)
	$(FINAL_CXX) $(CXXFLAGS) $(OBJ) -o $(TARGET) $(LDLIBS) $(MPI_LDLIBS)
	xgettext --keyword=_ --language=C --add-comments --sort-output -o po/mapgd.pot --copyright-holder="Matthew Ackerman" --package-name="mapgd" --package-version=$(MINOR) --msgid-bugs-address="matthew.s.ackerman@gmail.com"  $(addsuffix /*.cc, $(SUBDIR)) *.cc
	@echo "\n\nConsider typing 'make install DESTDIR=~/bin/' to install a local copy\n"

$(MPI_OBJ) : $(MPI_OBJ:.o=.cc)
	$(MPICXX) $(CXXFLAGS) -c $(@:.o=.cc) -o $@

DATA_SOURCE=data_types/*.cc

data.h : data_types/data.cc
	@echo "updating data_types"
	touch $(DATA_SOURCE)

interface.h : interface.cc
	@echo "updating commands"
	touch $(COMMAND_SOURCE)

debug : CXXFLAGS:=$(CXXFLAGS) $(DEBUG)
debug : all

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
	rm -f $(DOCDIR)/man/man1/*.gz
	for COM in $(COMMAND); do \
		if [ ! -f $(DOCDIR)/man/man1/mapgd-$$COM.1.gz ]; then \
			echo $(TARGET) $$COM \$$1 > mapgd-$$COM;  \
			chmod 777 mapgd-$$COM;  \
			help2man ./mapgd-$$COM > $(DOCDIR)/man/man1/mapgd-$$COM.1; \
			gzip $(DOCDIR)/man/man1/mapgd-$$COM.1; \
			rm mapgd-$$COM; \
		fi \
	done
	$(eval KEYS=`$(TARGET) keyinfo -a`) 
	echo -n $(KEYS)
	for KEY in $(KEYS); do \
		if [ ! -f $(DOCDIR)/man/man1/mapgd-$$TYPE.1.gz ]; then \
			echo $(TARGET) keyinfo $$KEY \$$1 > mapgd-$$KEY;  \
			chmod 777 mapgd-$$KEY;  \
			help2man ./mapgd-$$KEY > $(DOCDIR)/man/man1/mapgd-$$KEY.1; \
			gzip $(DOCDIR)/man/man1/mapgd-$$KEY.1; \
			rm mapgd-$$KEY; \
		fi \
	done
#.PHONY: clean
#	help2man

clean:
	rm -f objs/*.o
	rm -f objs/*.d
.PHONY: clean
