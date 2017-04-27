#ifndef _INDIVIDUAL_LIKELIHOOD_H_
#define _INDIVIDUAL_LIKELIHOOD_H_

#include <cmath>
#include <math.h>
#include <iomanip>      // std::setprecision
#include <cfloat>

#include "lnmultinomial.h"
#include "data_types/allele.h"
#include "newton-method-theta.h"
#include "models.h"

count_t init_params(Locus&, Allele&, const float_t&);										//!< Uses a huristic method for guessing good priors.
count_t maximize_analytical(Locus&, Allele&, models&, std::vector <float_t>&, const float_t&, const size_t&);	//!< Returns values assuming that all reads are identical.
count_t maximize_grid      (Locus&, Allele&, models&, std::vector <float_t>&, const float_t&, const size_t&);	//!< Uses a grid based search to maximize the likelihoood function. 
count_t maximize_newton    (Locus&, Allele&, models&, std::vector <float_t>&, const float_t&, const size_t&);	//!< Uses a newton-raphson method to find the maximum likelihoood priors. This function is still numerically unstable. 

count_t maximize_restricted_newton    (Locus&, Allele&, models&, std::vector <float_t>&, const float_t&, const size_t&);	//!< Uses a newton-raphson method to find the maximum likelihoood priors. This function is still numerically unstable. 
float_t loglikelihood(const Locus&, const Allele&, models&, const count_t&);	//!< Returns the loglikelihood of a set of quartets given the priors Allele

#endif
