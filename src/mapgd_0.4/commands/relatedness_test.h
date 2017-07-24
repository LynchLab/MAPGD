#ifndef _TEST_RELATEDNESS_H_
#define _TEST_RELATEDNESS_H_	

#include <stdio.h>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "interface.h" 
#include "map_file.h"
#include "genotype.h"
#include "genotype_pair.h"
#include "relatedness_data.h"
#include "relatedness.h"
//#include "newton-method-rel.h"

#ifndef NOOMP
#include <omp.h>
#endif

#ifndef NOGSL
#include <gsl/gsl_multimin.h>
#endif

/** \ingroup COMMANDS
  * @{*/
/// Estiamtes the relatedness of individuals with known genotypic probabilites.
int testRel(int argc, char *argv[]);
/** @}*/
#endif 
