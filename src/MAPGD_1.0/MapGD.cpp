#include <iostream>
#include <cstdio>
#include <cstring>
#include "Estimator.hpp"
#include "sam2pro.hpp"

using namespace std;

/** @brief Our main function.
  * Parses commandline arguments etc.
  * @returns 0 iff seccessful
**/
void version();
void usage();

int main (int argc, char* argv[]){

	const char *infile = "infile.txt";
	const char *outfile = "outfile.txt";

	if (argc<2) usage(); 
	if (std::strcmp(argv[1], "ep")== 0) {
		for (int optind=1; optind<argc; optind++){
			if (argv[optind][0]=='-'){
				char c=argv[optind][1];
				switch (c){
					case 0:break;
					case 'h': usage(); break;
					case 'v': version(); break;
					case 'i':
						infile=argv[optind+1];
						break;
					case 'o': 
						outfile=argv[optind+1];
						break;
					defalt :
						usage();
						break;
				}
			}
		}
		estimator(infile, outfile);
	}
	else if (std::strcmp(argv[1], "proview")== 0){
		sam2pro(argc-1, argv+1);
	}
	else usage();
};

void version(){
	const char str[]={
		"mapgd v1.0\n"
	};
	cout << str;
	exit(0);
};

void usage(){
	const char str[] ={
		"\nusage: mapgd <command> [OPTIONS].\n\n"
		"Command: ep\t\tEstimation allele frequency of pooled population data\n"
		"\t proview\tprint reads in the pro format\n"
		"\t -h\t\tPrint this.\n"
		"\t -v\t\tPrint version info.\n\n"
	};
	cout << str;
	exit(0);
};
