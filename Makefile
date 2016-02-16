SUBDIRS=src

$(SUBDIRS):
	$(MAKE) -C $@
.PHONY: $(SUBDIRS) 

nosql: 
	$(MAKE) -C $(SUBDIRS) nosql
.PHONY: nosql

all: $(SUBDIRS)
.PHONY: all

noomp: 
	$(MAKE) -C $(SUBDIRS) noomp
.PHONY: noomp 

docs: 
	$(MAKE) -C $(SUBDIRS) docs
.PHONY: docs

dist: 
	$(MAKE) -C $(SUBDIRS) dist
.PHONY: dist
	
test: all
	$(MAKE) -C $(SUBDIRS) test
.PHONY: test

install: all
	$(MAKE) -C $(SUBDIRS) all
.PHONY: install

clean:
	$(MAKE) -C $(SUBDIRS) clean
.PHONY: clean

push: 
	$(MAKE) -C $(SUBDIRS) all
	$(MAKE) -C $(SUBDIRS) test
	$(MAKE) -C $(SUBDIRS) clean
	git add -u
	git commit
	git push
.PHONY: push
