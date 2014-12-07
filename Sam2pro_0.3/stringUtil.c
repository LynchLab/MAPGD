/***** stringUtil.c **********************************************
 * Description: Collection of string handling functions.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Sun Jun  6 10:02:16 2004.
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
 ****************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "stringUtil.h"
#include "eprintf.h"

/* chomp: remove carriage return from string */
char *chomp(char *line){
    int i, l;

    l = strlen(line);

    for(i=0;i<l;i++){
      if(line[i] == '\n' || line[i] == '\r'){
	line[i] = '\0';
	break;
      }
    }
    return line;
}

/* fprintnf: print max of n characters of str onto fp; add ... if
 *   str was truncated
 */
void fprintnf(FILE *fp, char *str, int n){
  int i, l, m;
  
  l = strlen(str);
  m = n < l ? n : l;
  for(i=0;i<m;i++)
    fprintf(fp,"%c",str[i]);
  if(m < l)
    fprintf(fp,"...");
}

void strtolower(char *s, long l){
  long i;

  for(i=0;i<l;i++)
    s[i] = tolower(s[i]);
}

void strtoupper(char *s, long l){
  long i;

  for(i=0;i<l;i++)
    s[i] = toupper(s[i]);
}


/* strdup: make a duplicate of s */
char *strdup2(char *s){
  char *p;
  p = (char *) malloc(strlen(s)+1);    /* +1 for '\0' */
  if (p != NULL)
    p = strcpy(p,s);
  return p;
}

void split(char *line, char *splitC, char *splitArray[], int *arrayLen){
  char *line2, c;
  int i, j, start, end;

  start = end = 0;
  *arrayLen = 0;
  line2 = (char *) malloc(strlen(line)+1);    /* +1 for '\0' */
  if (line2 != NULL)
    line2 = strcpy(line2,line);
  i = 0;
  /* count number of fields */
  while((c = line2[end++]) != '\0'){
    if(c == *splitC){
      splitArray[*arrayLen] = (char *) malloc((end - start)*sizeof(char));
      i=0;
      for(j=start;j<end-1;j++){
	splitArray[*arrayLen][i++] = line2[j];
      }
      splitArray[*arrayLen][i] = '\0';
      start = end;
      (*arrayLen)++;
    }
  }
  i=0; 
  splitArray[*arrayLen] = (char *) malloc((end-start+1)*sizeof(char));
  for(j=start;j<end;j++){
    splitArray[*arrayLen][i++] = line2[j];
  }
  splitArray[*arrayLen][i] = '\0';
  (*arrayLen)++;
}

void replace(char *string, char original, char replacement){
  long i, l;
  l = strlen(string);
  for(i=0;i<l;i++){
    if(*string == original)
      *string = replacement;
    string++;
  }
}

/* remove leading and trailing non-alphanumerical characters */
char *cleanWordEdges(char *word){
  int l1, l2;
  char *word2 = NULL;

  /* remove leading non-alphanumericals */
  while(!isalnum(*word) && *word != '\0'){
    word++;
  }

  /* remove trailing non-alphanumericals */
  l1 = strlen(word);
  l2 = l1;
  while(l2 > 0 && !isalnum(word[l2-1])){
    l2--;
  }
  if(l2 < l1){
    word2 = (char *)emalloc((l2+1)*sizeof(char));
    word2 = strncpy(word2,word,l2);
    return word2; 
  }else{
    return word;
  }
  
}

/* remove all non-alphanumerical characters */
char *cleanWord(char *word){
  int i, j, l1, l2;
  char *word2 = NULL;

  /* get word length */
  l1 = strlen(word);
  
  /* count alphanumericals */
  l2=0;
  for(i=0;i<l1;i++)
    if(isalnum(word[i]))
      l2++;
  if(l2 < l1){  /* create cleaned string */
    word2 = (char *)emalloc((l2+1)*sizeof(char));
    j=0;
    for(i=0;i<l2;i++)
      if(isalnum(word[i]))
	 word2[j++] = word[i];
    word2[j] = '\0';
    return word2; 
  }else{
    return word;
  }
  
}

/* reverse: reverse stirng s in place */
void reverse(char *s){
  int c, i, j;
  for(i=0, j=strlen(s)-1; i<j; i++, j--){
    c = s[i];
    s[i] = s[j];
    s[j] = c;
  }
}

/* hash: form hash value for string s */
unsigned hash(char *s){
  unsigned hashval;
  
  for(hashval = 0; *s != '\0'; s++)
    hashval = *s +31 * hashval;

  return hashval % HASHSIZE;
}

/******************************** Word Counting with Binary Tree *******************************/
WordNode *addWord(WordNode *p, char *w){
  int cond;

  if(p == NULL){                       /* a new word has arrived */
    p = walloc();                      /* make a new node */
    p->word = strdup2(w);
    p->count = 1;
    p->left = p->right = NULL;
  }else if((cond = strcmp(p->word, w)) == 0)
    p->count++;       /* repeated word */
  else if(cond < 0)   /* less than into left subtree */
    p->left = addWord(p->left, w);
  else                /* greater than into right subtree */
    p->right = addWord(p->right, w);
  return p;
}

/* walloc: make WordNode */
WordNode *walloc(void){
  return (WordNode *) malloc(sizeof(WordNode));
}


/* treeprint: in-order print of tree p */
void treeprint(WordNode *p){
  if(p != NULL){
    treeprint(p->left);
    printf("%d %s\n",p->count,p->word);
    treeprint(p->right);
  }
}
/******************************** End of Word Counting with Binary Tree *******************************/
