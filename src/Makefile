SUBDIRS=stable
VER=1.0
NAME=mapgd

.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: subdirs

dist: 
	tar -czf $(NAME)-$(VER).tar.gz $(SUBDIRS)/*
	
test: all
	cd test && bash test.sh

install: all
	install -m 0755 ../bin/$(NAME) $(prefix)/bin
.PHONY: install

clean:
	cd $(SUBDIRS) && $(MAKE) $@
	rm -f $(NAME)-*.tar.gz
