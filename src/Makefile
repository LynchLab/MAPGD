SUBDIRS=mapgd_0.3.1
VER=0.3.1
NAME=mapgd
DEBIAN=0
UBUNTU=1

$(SUBDIRS):
	$(MAKE) -C $@
.PHONY: $(SUBDIRS)


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

install: all
	install -m 0755 ../bin/$(NAME) $(prefix)/bin
.PHONY: install

clean:
	cd $(SUBDIRS) && $(MAKE) $@
	rm -f $(NAME)-*.tar.gz
.PHONY: clean
