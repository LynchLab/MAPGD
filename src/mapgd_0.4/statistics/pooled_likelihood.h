#ifndef POOLEDLIKELIHOOD_H_
#define POOLEDLIKELIHOOD_H_

#include "lnmultinomial.h"
#include "../datatypes/allele.h"
#include <math.h>

void polymorphicmodel(allele const &, float_t *);
void monomorphicmodel(allele const &, float_t *);
void fixedmorphicmodel(allele const &, float_t *);

#endif
