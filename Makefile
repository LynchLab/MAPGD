SUBDIRS=src

.PHONY: clean
.PHONY: all
.PHONY: install
.PHONY: noomp

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: subdirs

noomp: all 
	?
dist: :
	tar -czf $(NAME)-$(VER).tar.gz $(SUBDIRS)/*
	
test: all
	cd test && bash test.sh

install: all
	install -m 0755 ./bin/$(NAME) $(prefix)/bin


clean:
	cd $(SUBDIRS) && $(MAKE) $@
	rm -f $(NAME)-*.tar.gz
