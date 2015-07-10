#ifndef INDLIKELIHOOD_H_
#define INDLIKELIHOOD_H_

#include <math.h>
#include "Likelihood.h"
#include "proFile.h"
#include <iomanip>      // std::setprecision
#include <cfloat>


float_t *mmModel(allele_stat *);
float_t *MmModel(allele_stat *);
float_t *MMModel(allele_stat *);

class models{
private:
lnmultinomial lnMM, lnMm, lnmm, lnF;	//The three multinomials we will use for probability calculations.
					//The '4' specifies the number of categories of the distribution.
					//Since these represent the distribution of the four nucleotides 
					//A, C, G and T, we use 4 categories.
public:
models(void);
float_t loglikelihood(site_t const &, allele_stat const &, count_t const &);
float_t lnP(count_t const *, allele_stat const &);
float_t genotypelikelihood(quartet_t const &, allele_stat const &, count_t const &);
};

count_t initparams(site_t &, allele_stat &, count_t const &, float_t const &, count_t const &);
count_t maximizegrid(site_t &, allele_stat &, models &, count_t const &, float_t const &, count_t const &);
float_t loglikelihood(site_t const &, allele_stat const &, models &, count_t const &);

/*@Breif, a function for calculating the jacobian of the log likelihood something something. */

#endif
