#include <iostream>
#include <cstdio>
#include "Estimator.hpp"

using namespace std;

/** @brief Our main function.
  * Parses commandline arguments etc.
  * @returns 0 iff seccessful
**/
void version();
void usage();

int main (int argc, char**argv){

	const char *infile = "infile.txt";
	const char *outfile = "outfile.txt";

	for (int optind=0; optind<argc; optind++){
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
};

void version(){
	const char str[]={
		"MapGD v1\n"
	};
	cout << str;
	exit(0);
};

void usage(){
	const char str[] ={
		"usage: MapGD -i [INFILE] -o [OUTFILE].\n"
		"\tINFILE needs to be a two column 'pro' file.\n"
		"Options:\n"
		"\t-h\tPrint this.\n"
		"\t-v\tPrint version info.\n"
	};
	cout << str;
	exit(0);
};
