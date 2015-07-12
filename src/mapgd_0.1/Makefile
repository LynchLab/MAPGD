CXXFLAGS= -g -std=c++11 -Wall -W -O3
prefix=/usr/local

mapdg: Estimator.o MapGD.o ReadFile.o sam2pro.o eprintf.o interface.o stringUtil.o tab.o
	$(CXX) $^ $(CXXFLAGS) -o ../../bin/mapgd

.PHONY: test
test:
	../scripts/RunMeLinux.sh ../test/test.pileup

.PHONY: debug
debug: CXXFLAGS=-g -std=c++11 -Wall -W -Og
debug: mapdg

    
install: mapdg
	install -m 0755 ../../bin/mapgd $(prefix)/bin
.PHONY: install

clean:
	rm -f *.o
	rm -f ../../bin/*
