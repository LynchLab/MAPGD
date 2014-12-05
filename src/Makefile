CXXFLAGS= -g -std=c++11 -Wall -W -O3

mapdg: Estimator.o MapGD.o 
	$(CXX) $^ $(CXXFLAGS) -o mapgd

.PHONY: test
test: CXXFLAGS+=-pg
test: clean mapdg
	./MapDG ../test/test.txt --

.PHONY: debug
debug: CXXFLAGS=-g -std=c++11 -Wall -W -Og 
debug: mapdg

clean:
	rm -f *.o mapgd
