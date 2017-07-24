/*  A command to compare two or more populations with pooled sequencing and 
 *  print Fst values.
 */
#ifndef _COMPAREPOOLED_H_
#define _COMPAREPOOLED_H_
#include "pro-file.h"
#include <iomanip>
#include "interface.h"
#include "compare-pooled.h"
#include "pooled_likelihood.h"
#include <iostream>
int comparePooled(int argc, char *argv[]);
#endif
