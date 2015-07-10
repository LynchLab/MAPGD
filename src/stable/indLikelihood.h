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
lnmultinomial *lnMM_, *lnMm_, *lnmm_, *lnF_;	//The three multinomials we will use for probability calculations.
					//The '4' specifies the number of categories of the distribution.
					//Since these represent the distribution of the four nucleotides 
					//A, C, G and T, we use 4 categories.
public:
models(void);
~models(void);
float_t loglikelihood(const site_t &, const allele_stat &, const count_t &);
float_t lnP(const count_t *, const allele_stat &);
float_t genotypelikelihood(const quartet_t &, const allele_stat &, const count_t &);
};

count_t initparams(site_t&, allele_stat&,const count_t&, const float_t&, const count_t &);
count_t maximizegrid(site_t&, allele_stat&, models&, std::vector <float_t>&, const count_t&, const float_t&, const count_t&);
float_t loglikelihood(const site_t&, const allele_stat&, models&, const count_t&);

#endif
