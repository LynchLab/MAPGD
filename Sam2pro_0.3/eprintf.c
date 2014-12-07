/***** eprintf.c **************************************************
 * Description: Collection of functions for error handling.
 * Reference: Kernighan, B. W. and Pike, R. (1999). The Practice
 *            of programming. Addision Wesley; chapter 4.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Fri Dec 17 11:16:34 2004.
 *
 * This file is part of kr.
 *
 *   kr is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   kr is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with kr; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *****************************************************************/
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include "eprintf.h"

/* efopen: open file and report if error */
FILE *efopen(char *fname, char *mode){
  FILE *fp;

  if((fp = fopen(fname, mode)) == NULL)
    eprintf("efopen(%s, %s) failed:",fname,mode);
  
  return fp;
}

/* eprintf: print error message and exit */
void eprintf(char *fmt, ...){
  va_list args;
  fflush(stdout);
  if(progname() != NULL)
    fprintf(stderr, "%s: ", progname());
  
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  
  if(fmt[0] != '\0' && fmt[strlen(fmt)-1] == ':')
    fprintf(stderr, " %s", strerror(errno));
  fprintf(stderr, "\n");
  exit(2); /* conventional value for failed execution */
}

/* estrdup: duplicate a string, report if error */
char *estrdup(char *s){
  char *t;
  
  t = (char *)malloc(strlen(s)+1);
  if(t == NULL)
    eprintf("estrdup(\"%.20s\") failed:", s);
  strcpy(t, s);
  return t;
}

/* emalloc: malloc and report if error */
void *emalloc(size_t n){
  void *p;
  
  p = malloc(n);
  if(p == NULL)
    eprintf("malloc of %u bytes failed:", n);
  return p;
}

/* erealloc: realloc and report if error */
void *erealloc(void *p, size_t n){

  p = realloc(p, n);
  if(p == NULL)
    eprintf("realloc of %u bytes failed:", n);
  return p;
}

static char *name = NULL; /* program name for messages */

/* setprogram: set stored name of program */
void setprogname2(char *str){
  name = estrdup(str);
}

/* progname: return stored name of program */
char *progname(void){
  return name;
}
