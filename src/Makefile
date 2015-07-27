SUBDIRS=mapgd_0.3.1
VER=0.3.1
NAME=mapgd


subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@
.PHONY: subdirs 


all: subdirs
.PHONY: all

noomp: 
	$(MAKE) -C $(SUBDIRS) noomp
.PHONY: noomp 

dist: 
	tar -czf $(NAME)-$(VER).tar.gz $(SUBDIRS)/*
	doxygen doxygen.conf
.PHONY: dist
	
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
