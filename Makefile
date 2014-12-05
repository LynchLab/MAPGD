SUBDIRS= src docs

.PHONY: subdirs $(SUBDIRS) code-docs clean

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

all: src

test: src
	cd test && sh testMapGD.sh

code-docs: subdirs
	$(MAKE) -C docs code-docs

clean:
	cd src && $(MAKE) $@
