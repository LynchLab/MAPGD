SUBDIRS=MAPGD_2.N 

.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: src

test: src
	cd test && sh testMapGD.sh
clean:
	cd src && $(MAKE) $@
