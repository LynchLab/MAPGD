SUBDIRS=src

$(SUBDIRS):
	$(MAKE) -C $@
.PHONY: $(SUBDIRS) 

all: $(SUBDIRS)
.PHONY: all

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
	echo $(HAVE_LIBM)
	echo $(HAVE_LIBGSLCBLAS)
.PHONY: echo
#define HAVE_LIBM 1
#define HAVE_LIBGSLCBLAS 1
#define HAVE_LIBGSL 1

	
test: all
	$(MAKE) -C $(SUBDIRS) test
.PHONY: test

install: all
	$(MAKE) -C $(SUBDIRS) all
.PHONY: install

clean:
	$(MAKE) -C $(SUBDIRS) clean
.PHONY: clean
	
increment:
	$(MAKE) -C $(SUBDIRS) increment 
.PHONY: clean

push: 
	$(MAKE) -C $(SUBDIRS) all 
	$(MAKE) -C $(SUBDIRS) test
	$(MAKE) -C $(SUBDIRS) clean
	$(MAKE) -C $(SUBDIRS) increment 
	git add -u
	git commit
	git push
.PHONY: push
