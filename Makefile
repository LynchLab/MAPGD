SUBDIRS=src

HAVE_GSL=no
HAVE_OMP=yes
HAVE_SQL=yes
HAVE_MAN=yes

CXX=g++
OPENMP_CXXFLAGS=-fopenmp
	
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
	@cat config.log
.PHONY: echo

	
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
	git checkout gh-pages
	./update.sh
	git checkout master
.PHONY: push
