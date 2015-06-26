#ifndef INDLIKELIHOOD_H_
#define INDLIKELIHOOD_H_

#include <math.h>
#include "Likelihood.h"
#include "proFile.h"
#include <iomanip>      // std::setprecision
#include <cfloat>

count_t initparams(site_t &, allele_stat &, count_t const &, float_t const &, count_t const &);
count_t maximizegrid(site_t &, allele_stat &, count_t const &, float_t const &, count_t const &);
float_t loglikelihood(site_t const &, allele_stat const &, count_t const &);
float_t *mmModel(allele_stat *);
float_t *MmModel(allele_stat *);
float_t *MMModel(allele_stat *);

#endif
