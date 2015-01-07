#include <iostream>

int main (int argc, char *argv[]){
	std::ostream *out=&std::cout;
	*out << "hi!\n";
};
