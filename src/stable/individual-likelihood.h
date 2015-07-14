#ifndef _INDIVIDUAL_LIKELIHOOD_H_
#define _INDIVIDUAL_LIKELIHOOD_H_

#include <math.h>

#include "lnmultinomial.h"
#include "allele_stat.h"

#include "pro-file.h"
#include <iomanip>      // std::setprecision
#include <cfloat>


struct args{
        bool verbose=false;
        bool quite=false;
	bool bayes=false;

        float_t minimum_error=0.001;
        float_t alpha=0.00;
        float_t maximum_gof=2.00;
        uint16_t maximum_excluded=96;
        uint16_t base_excluded=0;
        uint16_t minimum_coverage=0;
};

float_t *mmModel(allele_stat *);
float_t *MmModel(allele_stat *);
float_t *MMModel(allele_stat *);

/* \breif this class keeps the probabilities used in the calculation of linkelihoods.
*/
class models{
private:
	lnmultinomial *lnMM_, *lnMm_, *lnmm_, *lnF_;	//The three multinomials we will use for probability calculations.
							//The '4' specifies the number of categories of the distribution.
							//Since these represent the distribution of the four nucleotides 
							//A, C, G and T, we use 4 categories.
public:
	models(void);
	~models(void);
	float_t loglikelihood(const Locus &, const allele_stat &, const count_t &);
	float_t lnP(const count_t *, const allele_stat &);
	float_t genotypelikelihood(const quartet_t &, const allele_stat &, const count_t &);
};

count_t initparams(Locus&, allele_stat&, const args &);
count_t maximizegrid(Locus&, allele_stat&, models&, std::vector <float_t>&, const args&);
float_t loglikelihood(const Locus&, const allele_stat&, models&, const count_t&);

#endif
