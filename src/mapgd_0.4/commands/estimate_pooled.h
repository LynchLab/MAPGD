#ifndef _ESTIMATE_POOLED_H_
#define _ESTIMATE_POOLED_H_	

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "../mapgd.h"
#include "../statistics/pooled_likelihood.h"

int estimate_pooled(int, char **);
allele estimate(quartet_t);

#endif 
