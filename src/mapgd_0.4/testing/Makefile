#Testing with different standards.
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y

MODE=-O3
STD=-std=c++0x
WARNINGS=-Wall -Werror=return-type -Wno-write-strings
CXXFLAGS=$(STD) $(MODE) $(WARNINGS) -fopenmp 
DEBUG= -g 
CXX=g++


#>! TODO ORGANIZE THESE!!
OBJS=	../datatypes/row.o \
	../datatypes/formats.o \
	../datatypes/key.o \
	../datatypes/base_keys.o \
	key_test.o

key: $(OBJS)
	$(CXX) $^ $(CXXFLAGS) -o key_test

sql: sqlapi.o
	$(CXX) $^ $(CXXFLAGS) -o sql_test

clean:
	rm -f *.o
	rm -f key_test
	rm -f sql_test
