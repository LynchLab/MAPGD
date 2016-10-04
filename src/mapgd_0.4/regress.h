#ifndef _REGRESS_H_
#define _REGRESS_H_	

#include <stdio.h>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "interface.h" 
#include "map-file.h"
#include "genotype.h"
#include "genotype_pair.h"
#include "relatedness_data.h"

#ifndef NOOMP
#include <omp.h>
#endif

#ifndef NOGSL
#include <gsl/gsl_multimin.h>
#endif

/** \ingroup COMMANDS
  * @{*/
/// Estiamtes the population genetic components of phenotypic variation.
int regress(int argc, char *argv[]);
/** @}*/
#endif 
