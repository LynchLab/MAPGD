#ifndef _INDIVIDUAL_LIKELIHOOD_H_
#define _INDIVIDUAL_LIKELIHOOD_H_

#include <math.h>
#include <iomanip>      // std::setprecision
#include <cfloat>

#include "lnmultinomial.h"
#include "allele_stat.h"
#include "newton-method-theta.h"
#include "models.h"




count_t init_params(Locus&, allele_stat&, const float_t&);										//!< Uses a huristic method for guessing good priors.
count_t maximize_analytical(Locus&, allele_stat&, models&, std::vector <float_t>&, const float_t&, const size_t&);	//!< Returns values assuming that all reads are identical.
count_t maximize_grid      (Locus&, allele_stat&, models&, std::vector <float_t>&, const float_t&, const size_t&);	//!< Uses a grid based search to maximize the likelihoood function. 
count_t maximize_newton    (Locus&, allele_stat&, models&, std::vector <float_t>&, const float_t&, const size_t&);	//!< Uses a newton-raphson method to find the maximum likelihoood priors. This function is still numerically unstable. 
float_t loglikelihood(const Locus&, const allele_stat&, models&, const count_t&);	//!< Returns the loglikelihood of a set of quartets given the priors allele_stat

#endif
