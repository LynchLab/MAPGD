SUBDIRS= src src/Sam2pro_0.3

.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: src

test: src
	cd test && sh testMapGD.sh
clean:
	cd src && $(MAKE) $@
