#ifndef _RELATEDNESS_H_
#define _RELATEDNESS_H_	

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>

#include "interface.h" 
#include "map-file.h"
#include "genotype.h"
#include "relatedness_data.h"

void inc_f(Relatedness &, const Genotype_pair &, const size_t &);
void inc_theta(Relatedness &, const Genotype_pair &, const size_t &);
void inc_gamma(Relatedness &, const Genotype_pair &, const size_t &);
void inc_Delta(Relatedness &, const Genotype_pair &, const size_t &);

int relatedness(int, char **);

#endif 
