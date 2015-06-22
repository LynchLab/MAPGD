SUBDIRS=MAPGD_2.N 

.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: subdirs
	
test: all
	cd test && bash test.sh

install: all
	install -m 0755 ../bin/mapgd $(prefix)/bin
.PHONY: install

clean:
	cd src && $(MAKE) $@
