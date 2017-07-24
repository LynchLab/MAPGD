/*  A command to compare two or more populations with pooled sequencing and 
 *  print Fst values.
 */
#ifndef _COMPAREPOOLED_H_
#define _COMPAREPOOLED_H_

#include <iomanip>
#include <iostream>

#include "interface.h"
#include "compare-pooled.h"
#include "pooled_likelihood.h"

int comparePooled(int argc, char *argv[]);
#endif
