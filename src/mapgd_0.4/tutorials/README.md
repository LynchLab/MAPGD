First, compile avg\_coverage into *object* (https://en.wikipedia.org/wiki/Object_code) code (a .o file), which is just the machine code for ??? that has not been linked to the parts of mapgd ...
We will use the ... 

g++ -std=c++11 -O3 -fopenmp -c -o basic\_command.o basic\_command.cc

After this is done we need to *link* (https://en.wikipedia.org/wiki/Linker_%28computing%29) the object code 

g++ -o basic\_command basic\_command.o ../datatypes/\*.o
