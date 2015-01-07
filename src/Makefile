SUBDIRS= MAPGD_2.0 

.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: src

test: src
	cd test && sh testMapGD.sh
clean:
	cd src && $(MAKE) $@
