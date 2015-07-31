#include "convert.h"
#include <iostream>

void convert_file(const char *filename1, const char* filename2){
	std::cout << filename1 << "->" << filename2 << std::endl;
	streamtable(filename1, filename2);
};
