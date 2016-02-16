SUBDIRS=mapgd_0.4
VER=0.4
NAME=mapgd
DEBIAN=0
UBUNTU=1

$(SUBDIRS):
	$(MAKE) -C $@
.PHONY: $(SUBDIRS)

all: $(SUBDIRS)
.PHONY: all

nosql: 
	$(MAKE) -C $(SUBDIRS) noomp
.PHONY: nosql

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
	( cat doxygen.conf ; echo "PROJECT_NUMBER=$(VER)" ) | doxygen -
.PHONY: docs
	
test: all
	cd test && bash test.sh
.PHONY: test

clean: 
	-for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean); done
	rm -f $(NAME)-*.tar.gz
.PHONY: clean
