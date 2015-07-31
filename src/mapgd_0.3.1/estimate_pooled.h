#ifndef _ESTIMATE_POOLED_H_
#define _ESTIMATE_POOLED_H_	

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "interface.h" 
#include "pro_file.h" 
#include "pooled_likelihood.h"

int estimate_pooled(int, char **);
allele_stat estimate(quartet_t);

#endif 
