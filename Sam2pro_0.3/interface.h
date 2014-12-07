/***** interface.h ********************************
 * Description: Header file for user interface.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sun Jun 20 13:07:28 2004.
 * Licence: GNU General Public
 ***************************************************/ 
#define DEFAULT_M 4    /* minimum coverage */
#define DEFAULT_D "\t" /* tab delineator */


/* define argument container */
typedef struct args{
  char h;   /* help message? */
  char e;   /* error message? */
  char *d;  /* tab delineator */
  char **inputFiles;
  int m;    /* minimum coverage */
  int numInputFiles;
} Args;

Args *getArgs(int argc, char *argv[]);
void printUsage(char *version);
