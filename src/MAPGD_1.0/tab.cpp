/***** tab.c ******************************************************
 * Description: Implementation of tab library that handles tabular
 *              data. In contrast to the original description
 *              (c.f. ref.), this implementation can handle
 *              variable field separators.
 * Reference: Kernighan, B. W. and Pike, R. (1999). The Practice
 *            of programming. Addision Wesley; chapter 4.
 * Author: Bernhard Haubold, bernhard.haubold@fh-weihenstephan.de
 * File created on Wed Dec 22 09:56:03 2004.
 *****************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "tab.h"

enum {NOMEM = -2};                                /* out of memory signal */

/* constants & internal functions are declared 
 * static so they are visible only within the 
 * file that contains them.
 */
static char *line   = NULL;                       /* input chars */
static char *sline  = NULL;                       /* line copy used by split */
static int maxline  = 0;                          /* size of line[] and sline[] */
static char **field = NULL;                       /* field pointers */
static int maxfield = 0;                          /* size of field[] */
static int nfield   = 0;                          /* number of fields in field[] */
static char fieldsep[] = "\t";

static void reset(void);
static int endofline(FILE *fin, int c);
static int split(void);
static char *advquoted(char *p);
static char *getLine(FILE *fin);

/* setfieldsep: set field separator to an arbitrary character */
void tabSetFieldSep(char *c){
  fieldsep[0] = c[0];
}

/* tabGetPlainLine: get line without splitting in background */
char *tabGetPlainLine(FILE *fin){
  
  return getLine(fin);
  
}

static char *getLine(FILE *fin){
  int i, c;
  char *newl, *news;
  
  if(line == NULL){             /* allocate on first call */
    maxline = maxfield = 1;
    line = (char *) malloc(maxline);
    sline = (char *) malloc(maxline);
    field = (char **) malloc(maxfield*sizeof(field[0]));
    if(line == NULL || sline == NULL || field == NULL){
      reset();
      return NULL;               /* out of memory */
    } 
  }
  
  for(i=0; (c=getc(fin))!=EOF && !endofline(fin,c); i++){
    if(i >= maxline-1){          /* grow line */
      maxline *= 2;              /* double current size */
      newl = (char *) realloc(line, maxline);
      news = (char *) realloc(sline, maxline);
      if(newl == NULL || news == NULL){
	reset();
	return NULL;             /* out of memory */
      }
      line = newl;
      sline = news;
    }
    line[i] = c;
  }
  line[i] = '\0';

  return (c == EOF && i == 0) ? NULL : line;
}

/* tabGetLine: get one line, grow as needed
 * sample input:
 * "LU"\t86.25\t"11/4/98"\t"2:19PM"\t+4.0625
 */
char *tabGetLine(FILE *fin){
  char *newl;

  newl = getLine(fin);
  if(split() == NOMEM){
    reset();
    return NULL;                /* out of memory */
  }
  return newl;
}

static void reset(void){
  free(line);                   /* free(NULL) permitted by ANSI C */
  free(sline);
  free(field);
  line = NULL;
  sline = NULL;
  field = NULL;
  maxline = maxfield = nfield = 0;
}

/* endofline: check for and consume \r, \n, \r\n, or EOF */
static int endofline(FILE *fin, int c){
  int eol;

  eol = (c=='\r' || c=='\n');
  if(c == '\r'){
    c = getc(fin);
    if(c != '\n' && c != EOF)
      ungetc(c, fin);          /* read too far; put c back */
  }
  return eol;
}

/* split: split line into fields separated by fieldsep */
static int split(void){
  char *p, **newf;
  char *sepp;                 /* pointer to temporary separatpr character */
  int sepc;                   /* temporary separator character */
  
  nfield = 0;
  if(line[0] == '\0')
    return 0;
  strcpy(sline, line);
  p = sline;

  do{
    if(nfield >= maxfield){    
      maxfield *= 2;          /* double current size */
      newf = (char **) realloc(field, maxfield*sizeof(field[0]));
      if(newf == NULL)
	return NOMEM;
      field = newf;
    }
    if(*p == '"')
      sepp = advquoted(++p); /* skip initial quote */
    else
      sepp = p + strcspn(p, fieldsep);
    sepc = sepp[0];
    sepp[0] = '\0';         /* terminate field */
    field[nfield++] = p;
    p = sepp +1;
  }while(sepc == fieldsep[0]);
  
  return nfield;
}

/* advquoted: quoted field; return pointer to next separator */
static char *advquoted(char *p){
  int i, j, k;
  for(i=0,j=0; p[j] != '\0'; i++, j++){ 
    if(p[j] == '"' && p[++j] != '"'){    /* copy up to next separator or \0 */
      k = strcspn(p+j, fieldsep);
      i += k;
      j += k;
      break;
    }
  }
  p[i] = '\0';
  return p+j;
}

/* tabfield: return pointer to n-th field; counting starts at 0 */
char *tabField(int n){
  if(n<0 || n >= nfield)
    return NULL;
  return field[n];
}

/* tabnfield: return number of fields */
int tabNfield(void){
  return nfield;
}
     
