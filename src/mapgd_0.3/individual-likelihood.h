#ifndef _INDIVIDUAL_LIKELIHOOD_H_
#define _INDIVIDUAL_LIKELIHOOD_H_

#include <math.h>

#include "lnmultinomial.h"
#include "allele_stat.h"
#include "newton-method-theta.h"
#include "models.h"

#include "pro-file.h"
#include <iomanip>      // std::setprecision
#include <cfloat>



count_t init_params(Locus&, allele_stat&,const count_t&, const float_t&, const count_t &);
count_t maximize_analytical(Locus&, allele_stat&, models&, std::vector <float_t>&, const count_t&, const float_t&, const count_t&);
count_t maximize_grid      (Locus&, allele_stat&, models&, std::vector <float_t>&, const count_t&, const float_t&, const count_t&);
count_t maximize_newton    (Locus&, allele_stat&, models&, std::vector <float_t>&, const count_t&, const float_t&, const count_t&);
float_t loglikelihood(const Locus&, const allele_stat&, models&, const count_t&);

#endif
