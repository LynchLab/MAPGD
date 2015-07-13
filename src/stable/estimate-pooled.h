#ifndef _ESTIMATE_POOLED_H_
#define _ESTIMATE_POOLED_H_	

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "interface.h" 
#include "pro-file.h" 
#include "pooled-likelihood.h"

int estimatePooled(int, char **);
allele_stat estimate(quartet_t);

#endif 
