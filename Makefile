SUBDIRS=src

PACKAGE_TARNAME = mapgd

prefix      = /usr/local
exec_prefix = ${prefix}
bindir      = ${exec_prefix}/bin
includedir  = ${prefix}/include
libdir      = ${exec_prefix}/lib
libexecdir  = ${exec_prefix}/libexec
datarootdir = ${prefix}/share
mandir      = ${datarootdir}/man
man1dir     = ${datarootdir}/man/man1
man5dir     = ${datarootdir}/man/man5
DOCPATH     = ${datarootdir}/doc/${PACKAGE_TARNAME}

HAVE_GSL=yes
HAVE_OMP=yes
HAVE_MPI=
HAVE_SQL=yes
HAVE_HTS=yes
HAVE_MAN=yes
USE_NLS=yes
HAVE_LZ4=@HAVE_LZ4@
PACKAGE_VERSION=0.4.38 
PROGRAM_VERSION=$(shell cat ./src/mapgd_0.4/VERSION | cut -d '-' -f 1)

export OMPI_CXX=g++

MPICXX=g++ -std=c++11
CXX=g++
OPENMP_CXXFLAGS=-fopenmp
CXXFLAGS=-g -O2
LDFLAGS=

LDLIBS += 

#SAMTOOLS FLAGS

##HTSDIR = 
##HTSLIB = $(HTSDIR)/libhts.a
##HTSLIB_LIB = $(HTSLIB)
##BGZIP = $(HTSDIR)/bgzip

# HTSDIR is not always defined in ax_with_htslib
#HTSDIR := $(realpath )
HTSLIB_CPPFLAGS = 
HTSLIB_LDFLAGS =  
HTSLIB_LIB = -lhts

#LZ4 FLAGS

LZ4_LIB = -L  -l llz4
LZ4_INCLUDE = -I 

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


#test needs to prompt and remove untracked files.
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
	PROGRAM_VERSION=$(shell cat ./src/mapgd_0.4/VERSION)
	git tag $(PROGRAM_VERSION)
	git add -u
	git commit
	git push 
	git push origin $(PROGRAM_VERSION)
	rm -rf autom4te.cache/
	rm -rf aclocal.m4
	rm -rf config.log
	rm -rf config.status
	git checkout gh-pages
	./update.sh
	git checkout master
.PHONY: push
