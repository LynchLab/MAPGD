STATGENDIR=../libStatGen_1.0.8/
LDLIBS= -lStatGen -lz 
LDFLAGS= -L $(STATGENDIR) 
 
CXXFLAGS= -I -g -std=c++0x -Wall -W $(LDLIBS) $(LDFLAGS) -O3 -fstack-protector -fopenmp -pthread

dgrid: dgrid.o compute.o genotype.o readtable.o convert.o
	$(CXX) $^ $(CXXFLAGS) -o dgrid

.PHONY: prof
prof: CXXFLAGS+=-pg
prof: clean dgrid
	./dgrid ../test/test.txt 1 2

.PHONY: debug
debug: CXXFLAGS= -I $(STATGENDIR)/include -g -std=gnu++0x -Wall -W $(LDLIBS) $(LDFLAGS) -Og -fstack-protector-all
debug: dgrid

clean:
	rm -f *.o dgrid
