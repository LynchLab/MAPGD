#include <iostream>
#include <cstdio>
#include <getopt.h>

#include "compute.hpp"
#include "convert.hpp"

#define THREADS 32

using namespace std;

/** @brief Our main function.
  * Parses commandline arguments, opens the file and starts computation.
  * @returns 0 iff seccessful
**/
void version();
void usage();

int main (int argc, char**argv){
        static struct option long_options[] = {
                {"version", no_argument, NULL, 'v'},
                {"help", no_argument, NULL, 'h'},
                {"A", required_argument, NULL, 'A'},
                {"B", required_argument, NULL, 'B'},
                {"min", required_argument, NULL, 'm'},
                {"max", required_argument, NULL, 'M'},
                {"user", required_argument, NULL, 'U'},
                {"subsample", required_argument, NULL, 's'},
                {"convert", required_argument, NULL, 'c'},
                {NULL,0,NULL,0}
        };

	size_t a=0; 
	size_t b=1;

	size_t m=0; 
	size_t M=1000;

	const char *filename = nullptr;
	const char *filename2 = nullptr;
	const char *pFile = nullptr;
	bool cvt=false;
	size_t subsample=0;

	while (1){
		int option_index = 0;
		int c = getopt_long(argc, argv, "h:A:B:s:m:M:v:U:c", long_options, &option_index);

		if( c == -1){
                        break;
                }

		switch (c){
			case 0:break;
			case 'h': usage(); break;
			case 'c': 
				cvt=true;
				filename2=argv[optind];
				break;
			case 'A':
				a = stoull(optarg);
				break;
			case 'B':
				b = stoull(optarg);
				break;
                        case 'm':
                                m = stoull(optarg);
                                break;
                        case 'M':
                                M = stoull(optarg);
                                break;
                        case 'U':
                                pFile = optarg;
                                break;
                        case 's':
                                subsample = stoull(optarg);
                                break;
			case '?': 
			defalt :
				usage();
				break;
		}
	}
	if (optind==argc){
		cerr << "Missing filename. Aborting." << endl;
		exit(1);
	}
	else filename=argv[optind];
	if (cvt) convert_file(filename, filename2);
	else compute(filename, a, b, subsample);
};

void version(){
	const char str[]={
		"dgrid v e^(i*p)"
	};
	cout << str;
	exit(0);
};

void usage(){
	const char str[] ={
		"usage: dgrid [OPTIONS] FILE.\n"
		"\tFILE needs to be a parsed takahiro file.\n"
		"Options:\n"
		"\t-A\tThe Ath individual in the population.\n"
		"\t-B\tThe Bth individual in the population.\n"
		"\t-h\tPrint this.\n"
		"\t-v\tPrint version info (don't do that!).\n"
	};
	cout << str;
	exit(0);
};
