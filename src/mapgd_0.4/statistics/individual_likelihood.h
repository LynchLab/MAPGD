#ifndef _INDIVIDUAL_LIKELIHOOD_H_
#define _INDIVIDUAL_LIKELIHOOD_H_

#include <math.h>
#include <iomanip>      // std::setprecision
#include <cfloat>

#include "../datatypes/allele.h"
#include "lnmultinomial.h"
#include "newton_method_theta.h"
#include "models.h"
#include "../data.h"

count_t init_params(locus&, allele&, const float_t&);										//!< Uses a huristic method for guessing good priors.
count_t maximize_analytical(locus&, allele&, models&, std::vector <float_t>&, const float_t&, const size_t&);	//!< Returns values assuming that all reads are identical.
count_t maximize_grid      (locus&, allele&, models&, std::vector <float_t>&, const float_t&, const size_t&);	//!< Uses a grid based search to maximize the likelihoood function. 
count_t maximize_newton    (locus&, allele&, models&, std::vector <float_t>&, const float_t&, const size_t&);	//!< Uses a newton-raphson method to find the maximum likelihoood priors. This function is still numerically unstable. 
float_t loglikelihood(const locus&, const allele&, models&, const count_t&);	//!< Returns the loglikelihood of a set of quartets given the priors allele
real_t efc (const locus &site);							//!< returns the number of effective chromosomes at a locus

#endif
