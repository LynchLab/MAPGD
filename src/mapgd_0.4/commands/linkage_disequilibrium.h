#ifndef _POPLD_H_
#define _POPLD_H_	

#include "typedef.h"
#include "interface.h"
#include "likelihood.h"
#include "genotype.h"
#include "locus.h"
#include "map_file.h"
#include "newton-method-ld.h"
#include "datatypes.h"

#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <queue>
#include <algorithm>
#include <limits.h>
#include <iterator>     // std::distance

#ifndef NOOMP
#include <omp.h>
#endif 

int linkage_disequilibrium(int, char **);

#endif 
