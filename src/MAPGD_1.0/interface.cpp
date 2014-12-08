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
#include <unistd.h>
#include "interface.hpp"
#include "eprintf.hpp"

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
	printf("ERROR: please specify 5 or 6 column output\n");
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
      printf("# unknown argument: %c\n",c);
      args->e = 1;
      return args;
    }
    c = getopt(argc, argv, optString);
  }
  args->inputFiles = argv + optind;
  args->numInputFiles = argc - optind;
  if(args->c == 0){
    printf("ERROR: please specify the desired number of output columns using the \ty{-c} option\n");
    args->e = 1;
  }
  return args;
}


void printUsage(char *version){
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

