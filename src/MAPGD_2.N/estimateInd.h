#ifndef _ESTIMATOR_IND_H_
#define _ESTIMATOR_IND_H_	

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "interface.h" 
#include "proFile.h" 
#include "indLikelihood.h"
#include "Likelihood.h"

int estimateInd(int, char **);
allele_stat_t estimate(quartet_t);
#endif 
