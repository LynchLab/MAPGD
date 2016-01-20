SUBDIRS=pyer
VER=0.3.1
NAME=mapgd

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: $(SUBDIRS)
.PHONY: all

clean:
	-for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean ); done
.PHONY: clean
