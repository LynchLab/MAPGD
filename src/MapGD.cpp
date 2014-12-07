#include <iostream>
#include <cstdio>
//#include <getopt.h>

#include "Estimator.hpp"

using namespace std;

/** @brief Our main function.
  * Parses commandline arguments etc.
  * @returns 0 iff seccessful
**/
void version();
void usage();

int main (int argc, char**argv){
	/*
        static struct option long_options[] = {
                {"version", no_argument, NULL, 'v'},
                {"help", no_argument, NULL, 'h'},
                {"outfile", required_argument, NULL, 'o'},
                {"infile", required_argument, NULL, 'i'},
                {NULL,0,NULL,0}
        };*/

	const char *infile = "infile.txt";
	const char *outfile = "outfile.txt";

	/*while (1){
		int option_index = 0;
		int c = getopt_long(argc, argv, "hv::oi", long_options, &option_index);

		if( c == -1){
                        break;
                }
*/
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
