VER=0.4
SUBDIRS=mapgd_$(VER)
NAME=mapgd
DEBIAN=0
UBUNTU=1

export VER 

$(SUBDIRS):
	$(MAKE) -C $@
.PHONY: $(SUBDIRS)

all: $(SUBDIRS)
	echo $(NOOMP) $(NOSQL)
.PHONY: all

nosql: 
	$(MAKE) -C $(SUBDIRS) nosql
.PHONY: nosql

noomp: 
	$(MAKE) -C $(SUBDIRS) noomp
.PHONY: noomp

increment: 
	$(MAKE) -C $(SUBDIRS) increment
.PHONY: increment

dist: docs
	tar -czf $(NAME)-$(VER).tar.gz $(SUBDIRS)/*
	bzr dh-make $(NAME) $(VER) $(NAME)-$(VER).tar.gz
	rm -rf debian/*ex debian/*EX
	sed $(VER)-$(DEBIAN)ubuntu$(UBUNTU)	
.PHONY: dist

docs: 
	$(MAKE) -C $(SUBDIRS) docs
	( cat doxygen.conf ; echo "PROJECT_NUMBER=$(VER)" ) | doxygen -
.PHONY: docs
	
test: all
	cd test && bash test.sh
.PHONY: test

clean: 
	$(MAKE) -C $(SUBDIRS) clean
	rm -f $(NAME)-*.tar.gz
.PHONY: clean
