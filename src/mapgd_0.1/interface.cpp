/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:12:10 2004.
 *****************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <iostream>

#include "interface.hpp"
#include "eprintf.hpp"

extern char *optarg;
extern int optind, opterr, optopt;

/* @Breif : A reimplementation of the unix getopt command for compatibility.
*  It probably doesn't follow all the standards yet and will no doubt break things.
*/

char getopt(int argc, char *argv[], char*optString){
	char *opt=optString;

	if (optind<argc){
		if (argv[optind][0]=='-'){
			optopt=argv[optind][1];
			while(*opt!=0){
				if (optopt==*opt){
					if ( *(opt+1)==':') if (optind<argc-1) optarg=argv[optind+1];
					optind+=2;
					return *opt;
				};
				opt++;
			};
			optind+=1;
			return '?';
		};
		return -1;
	};
	return -1;
};


Args *args;

Args *getArgs(int argc, char *argv[]){
  char c;
  char *optString = "hm:d:c:";

  args = (Args *)emalloc(sizeof(Args));
  args->m = DEFAULT_M;
  args->d = DEFAULT_D;
  args->c = 0;
  args->h = 0;
  args->e = 0;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
    case 'm':                           /* minimum coverage */
      args->m = atoi(optarg);
      break;
    case 'c':                           /* number of columns in output */
      args->c = atoi(optarg);
      if(args->c != 5 && args->c != 6){
	std::cerr << "ERROR: please specify 5 or 6 column output\n";
	args->e = 1;
      }
      break;
    case 'd':
      args->d = optarg;
      break;
    case '?':                           /* fall-through is intentional */
    case 'h':                           /* print help */
      args->h = 1;
      break;
    default:
      std::cerr << "# unknown argument:" << c << std::endl;
      args->e = 1;
      return args;
    }
    c = getopt(argc, argv, optString);
  }
  args->inputFiles = argv + optind;
  args->numInputFiles = argc - optind;
  if(args->c == 0){
    std::cerr << "ERROR: please specify the desired number of output columns using the \ty{-c} option" << std::endl;
    args->e = 1;
  }
  return args;
}


void printUsage(char const *version){
  printf("sam2pro version %s written by Bernhard Haubold\n", version);
  printf("purpose: convert sam output to profiles\n");
  printf("usage: sam2pro [inputFile(s)]\n");
  printf("options:\n");
  printf("\t-c <NUM> number of columns in output: 5|6]\n");
  printf("\t[-m <NUM> minimum coverage; default: NUM=%d]\n",DEFAULT_M);
  printf("\t[-d <CHAR> tab delineator; default: CHAR=%s]\n",DEFAULT_D);
  printf("\t[-h print this help message and exit]\n");
  exit(0);
}

