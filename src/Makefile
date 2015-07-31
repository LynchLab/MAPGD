SUBDIRS=mapgd_0.3.1 intefaces 
VER=0.3.1
NAME=mapgd
DEBIAN=0
UBUNTU=1


.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: $(SUBDIRS)
.PHONY: all

noomp: 
	$(MAKE) -C $(SUBDIRS) noomp
.PHONY: noomp 

dist: docs 
	tar -czf $(NAME)-$(VER).tar.gz $(SUBDIRS)/*
	bzr dh-make $(NAME) $(VER) $(NAME)-$(VER).tar.gz
	rm -rf debian/*ex debian/*EX
	sed $(VER)-$(DEBIAN)ubuntu$(UBUNTU)	
.PHONY: dist

docs: 
	doxygen doxygen.conf
.PHONY: docs
	
test: all
	cd test && bash test.sh
.PHONY: test

clean: 
	-for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean ); done
	rm -f $(NAME)-*.tar.gz
.PHONY: clean
