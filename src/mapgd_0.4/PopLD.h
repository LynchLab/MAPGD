#ifndef _POPLD_H_
#define _POPLD_H_	

#include "typedef.h"
#include "interface.h"
#include "likelihood.h"
#include "genotype.h"
#include "locus.h"
#include "map-file.h"
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

#ifndef NOOMP
#include <omp.h>
#endif 

int PopLD(int, char **);

#endif 
