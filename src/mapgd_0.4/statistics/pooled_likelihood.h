#ifndef POOLEDLIKELIHOOD_H_
#define POOLEDLIKELIHOOD_H_

#include "lnmultinomial.h"
#include "../datatypes/allele.h"
#include <math.h>

void polymorphicmodel(allele_t const &, float_t *);
void monomorphicmodel(allele_t const &, float_t *);
void fixedmorphicmodel(allele_t const &, float_t *);

#endif
