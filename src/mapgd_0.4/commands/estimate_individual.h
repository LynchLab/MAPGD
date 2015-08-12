#ifndef _ESTIMATE_INDIVIDUAL_H_
#define _ESTIMATE_INDIVIDUAL_H_	

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>

#include <iostream>
#include <iomanip>
#include <ctime>

#include <algorithm>
#include <functional>

#include "../interface.h" 
#include "../statistics/individual_likelihood.h"
#include "../statistics/likelihood.h"
#include "../pipe/pipe.h"
#include "../data.h"

///The estimate_individual command. Commands need to function like main. 
int estimate_individual(int, char **);

row estimate (locus &site, models &model, const count_t &MIN, const float_t &EMLMIN, const float_t &MINGOF, const size_t &MAXPITCH);

#endif 
