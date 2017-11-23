#ifndef _REML_H_
#define _REML_H_	

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
#include "relatedness.h"

#ifndef NOOMP
#include <omp.h>
#endif

#ifdef MPI
#include <ciso646>
#include <mpi.h>
#endif

int reml(int, char **);

#endif 
