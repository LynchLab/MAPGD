#ifndef _RELATEDNESS_H_
#define _RELATEDNESS_H_	

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include <gsl/gsl_multimin.h>

#include "interface.h" 
#include "map-file.h"
#include "genotype.h"
#include "genotype_pair.h"
#include "relatedness_data.h"

#include <omp.h>

void inc_f(Relatedness &, const Genotype_pair &, const size_t &);
void inc_theta(Relatedness &, const Genotype_pair &, const size_t &);
void inc_gamma(Relatedness &, const Genotype_pair &, const size_t &);
void inc_Delta(Relatedness &, const Genotype_pair &, const size_t &);

int estimateRel(int argc, char *argv[]);

#endif 
