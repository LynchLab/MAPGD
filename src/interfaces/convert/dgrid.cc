/* This is a snipit of an old program used to due relatedness by itself. It's just around for converting a file will be
 * ditched soon. */

#include <iostream>
#include <cstdio>
#include <getopt.h>

#include "convert.h"

int main (int argc, char**argv){
	if (argc!=3){
		cerr << "Usage FILE1 FILE2\nconverts FILE1 to FILE2, which is readable by relatedness.py" << endl;
		exit(1);
	}
	convert_file(argv[1], argv[2]);
}
