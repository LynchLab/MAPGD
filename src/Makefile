SUBDIRS=mapgd_0.3.1
VER=0.3.1
NAME=mapgd

.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: subdirs

noomp: subdirs noomp

dist: 
	tar -czf $(NAME)-$(VER).tar.gz $(SUBDIRS)/*
	doxygen doxygen.conf
	
test: all
	cd test && bash test.sh

install: all
	install -m 0755 ../bin/$(NAME) $(prefix)/bin
.PHONY: install

clean:
	cd $(SUBDIRS) && $(MAKE) $@
	rm -f $(NAME)-*.tar.gz
