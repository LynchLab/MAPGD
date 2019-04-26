#ifndef _RELATEDNESS_H_
#define _RELATEDNESS_H_	

#include <stdio.h>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
//#include <random>

#include "interface.h" 
#include "map_file.h"
#include "genotype.h"
#include "genotype_pair.h"
#include "relatedness_data.h"
#include "newton-method-rel.h"
#include "constants.h"

#ifndef NOOMP
#include <omp.h>
#endif

#ifndef NOGSL
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#endif

/** \ingroup COMMANDS
  * @{*/
/// Estiamtes the relatedness of individuals with known genotypic probabilites.
int estimateRel(int argc, char *argv[]);
/** @}*/

size_t
freqtoi(float_t in);

//DONE Moved to in memory
std::map <Genotype_pair_tuple, size_t> 
hash_genotypes (const std::stringstream &, const size_t &, const size_t &, const bool &, const bool &);

//DONE Moved to in memory
std::map <Genotype_pair_tuple, size_t> 
downsample_genotypes (const std::stringstream &, const size_t &, const size_t &, const bool &, const bool&);

/*Does a regression of allele frequency of the samples on the population allele frequency*/
void 
set_e(Relatedness &relatedness, std::map <Genotype_pair_tuple, size_t> &hashed_genotypes);

/*Guess starting values of relatedness for the maximization procedure*/
void 
gestimate(Relatedness &relatedness, std::map <Genotype_pair_tuple, size_t> &counts);

void 
inc_f(Relatedness &rel, const Genotype_pair &pair, const size_t &count);

void 
inc_theta(Relatedness &rel, const Genotype_pair &pair, const size_t &count);

void 
inc_gamma(Relatedness &rel, const Genotype_pair &pair, const size_t &count);

void 
inc_Delta (Relatedness &rel, const Genotype_pair &pair, const size_t &count);

//TODO FIX THIS TERRIBLE HACK JOB YOU SCHMUCK!

float_t
get_ll (const Relatedness &rel, const Genotype_pair &pair, const float_t count);

#ifndef NOGSL
double
rel_ll (const gsl_vector *v, void *void_hashed_genotypes_p);
#endif 

/*Maximizes the relatedness*/
void 
maximize(Relatedness &rel, std::map <Genotype_pair_tuple, size_t> &hashed_genotypes);

int
newton(Relatedness &rel, std::map <Genotype_pair_tuple, size_t> &hashed_genotypes);

void
get_llr(Relatedness &rel, std::map <Genotype_pair_tuple, size_t> hashed_genotypes);

#endif 

