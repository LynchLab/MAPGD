/*  A command to compare two or more populations with pooled sequencing and 
 *  print Fst values.
 */
#ifndef _COMPAREPOOLED_H_
#define _COMPAREPOOLED_H_

#include "../mapgd.h"
#include "../statistics/pooled_likelihood.h"

int compare_pooled(int argc, char *argv[]);
#endif
