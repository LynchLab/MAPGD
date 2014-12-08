/***** interface.c ************************************************
 * Description: Routine for gathering arguments from the command
 *              line.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:12:10 2004.
 *****************************************************************/

/*Minor modifications by Matt Ackerman to run under windows*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#include "interface.h"
#include "eprintf.h"

Args *args;

Args *getArgs(int argc, char *argv[]){
  char c;
  char *optString = "hm:d:";

  args = (Args *)emalloc(sizeof(Args));
  args->m = DEFAULT_M;
  args->d = DEFAULT_D;
  args->h = 0;
  args->e = 0;

  c = getopt(argc, argv, optString);
  while(c != -1){
    switch(c){
    case 'm':                           /* minimum coverage */
      args->m = atoi(optarg);
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
  return args;
}


void printUsage(char *version){
  printf("sam2pro version %s written by Bernhard Haubold\n", version);
  printf("purpose: convert sam output to profiles\n");
  printf("usage: sam2pro [inputFile(s)]\n");
  printf("options:\n");
  printf("\t[-m <NUM> minimum coverage; default: NUM=%d]\n",DEFAULT_M);
  printf("\t[-d <CHAR> tab delineator; default: CHAR=%s]\n",DEFAULT_D);
  printf("\t[-h print this help message and exit]\n");
  exit(0);
}

