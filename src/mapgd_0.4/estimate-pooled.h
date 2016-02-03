#ifndef _ESTIMATE_POOLED_H_
#define _ESTIMATE_POOLED_H_	

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "interface.h" 
#include "map-file.h"
#include "pooled-likelihood.h"
#include "pooled_data.h"
#include "locus.h"
#include "allele_stat.h"

int estimatePooled(int, char **);
allele_stat estimate(quartet_t);

#endif 
