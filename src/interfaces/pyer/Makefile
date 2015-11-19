VER=0.3.1
NAME=mapgd

all:
	python setup.py build
	rm -rf ../../../lib/rml.so
	ln -sf `pwd`/`find build -name *.so` ../../../lib/ 
.PHONY: all

clean:
	rm -rf build 
	rm -rf ../../../lib/rml.so
.PHONY: clean
