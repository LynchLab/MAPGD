SUBDIRS=src

HAVE_GSL=yes
HAVE_OMP=yes
HAVE_MPI=yes
HAVE_SQL=yes
HAVE_HTS=yes
HAVE_MAN=yes

MPICXX=g++ -std=c++11 -std=c++11
CXX=g++ -std=c++11
OPENMP_CXXFLAGS=-fopenmp
CXXFLAGS=-O4

#SAMTOOLS FLAGS

HTSDIR = ../htslib
#include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
HTSLIB_LIB = $(HTSLIB)
BGZIP = $(HTSDIR)/bgzip
HTSLIB_CPPFLAGS = -I ../../../htslib
HTSLIB_LDFLAGS = -L ../../../htslib
#HTSLIB_LIB = -lhts
	
export

all:  
all: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@
.PHONY: $(SUBDIRS) 

nosql: NOSQL = 1
nosql: all
.PHONY: nosql


noomp: NOOMP = 1  
noomp: all
.PHONY: noomp 

docs: 
	$(MAKE) -C $(SUBDIRS) docs
.PHONY: docs

dist: 
	$(MAKE) -C $(SUBDIRS) dist
.PHONY: dist

echo:
	@echo HAVE_GSL $(HAVE_GSL)
	@echo HAVE_OMP $(HAVE_OMP)
	@echo HAVE_SQL $(HAVE_SQL)
	@echo HAVE_HTS $(HAVE_HTS)
	@cat config.log
.PHONY: echo

debug: 
	$(MAKE) -C $(SUBDIRS) debug
.PHONY: dist
	
test: all
	$(MAKE) -C $(SUBDIRS) test
.PHONY: test

install: all
	$(MAKE) -C $(SUBDIRS) install
.PHONY: install

clean:
	rm -f config.log  
	rm -f config.status
	$(MAKE) -C $(SUBDIRS) clean
.PHONY: clean
	
increment:
	$(MAKE) -C $(SUBDIRS) increment 
.PHONY: clean

push: 
	$(MAKE) -C $(SUBDIRS) all 
	$(MAKE) -C $(SUBDIRS) test
	$(MAKE) -C $(SUBDIRS) docs
	$(MAKE) -C $(SUBDIRS) clean
	rm -f ./bin/mapgd
	CXX=i686-w64-mingw32-c++
	TARGET=mapgd-win64.exe
	$(MAKE) -C $(SUBDIRS) all 
	$(MAKE) -C $(SUBDIRS) clean
	CXX=x86_64-w64-mingw32-g++
	TARGET=mapgd-win32.exe
	$(MAKE) -C $(SUBDIRS) all 
	$(MAKE) -C $(SUBDIRS) clean
	rm -f ./bin/mapgd
	cp Makefile.empty Makefile
	$(MAKE) -C $(SUBDIRS) increment 
	git add -u
	git commit
	git push
	rm -rf autom4te.cache/
	git checkout gh-pages
	./update.sh
	git checkout master
.PHONY: push
