#ifndef _ESTIMATOR_IND_H_
#define _ESTIMATOR_IND_H_	

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <ctime>
#include <chrono>

#include "interface.h" 
#include "proFile.h" 
#include "indLikelihood.h"
#include "Likelihood.h"
//#include "pgdFile.h"

int estimateInd(int, char **);
allele_stat estimate(quartet_t, count_t, count_t, count_t, count_t);
#endif 
