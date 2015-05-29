#ifndef _ESTIMATOR_HPP_
#define _ESTIMATOR_HPP_	

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "interface.h" 
#include "proFile.h" 
#include "pooledLikelihood.h"

int estimatePooled(int, char **);

allele_stat_t estimate(quartet_t);
#endif 
