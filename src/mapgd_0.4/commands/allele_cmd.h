#ifndef _ALLELE_CMD_H_
#define _ALLELE_CMD_H_	

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <ctime>

#include <algorithm>
#include <functional>

#include "interface.h" 
#include "map_file.h"
#include "individual_likelihood.h"
#include "likelihood.h"
#include "genotype.h"
#include "sample_gof.h"
#include "locus.h"

#ifndef NOOMP
#include <omp.h>
#endif

#ifdef MPI
#include <ciso646>
#include <mpi.h>
#endif

int allele_cmd(int, char **);
Allele estimate (Locus &site, models &model, std::vector<float_t> &gofs, const count_t &MIN, const float_t &EMLMIN, const float_t &MINGOF, const size_t &MAXPITCH);

#endif 
