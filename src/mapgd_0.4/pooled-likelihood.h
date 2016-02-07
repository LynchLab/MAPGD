#ifndef POOLEDLIKELIHOOD_H_
#define POOLEDLIKELIHOOD_H_

#include "lnmultinomial.h"
#include "data_types/allele.h"
#include <math.h>

void polymorphicmodel(Allele const &, float_t *);
void monomorphicmodel(Allele const &, float_t *);
void fixedmorphicmodel(Allele const &, float_t *);

#endif
