#ifndef _ESTIMATE_POOLED_H_
#define _ESTIMATE_POOLED_H_	

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "interface.h" 
#include "map_file.h"
#include "pooled_likelihood.h"
#include "pooled_data.h"
#include "locus.h"
#include "data_types/allele.h"

/** 
  * \ingroup COMMANDS
  * @{
  */
///Estimates allele fequency statistics from pooled data
int estimatePooled(int, char **);
/** @}*/

Allele estimate(quartet_t);

#endif 
