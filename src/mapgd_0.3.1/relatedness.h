#ifndef _RELATEDNESS_H_
#define _RELATEDNESS_H_

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
#include "pro-file.h" 
#include "individual-likelihood.h"
#include "likelihood.h"
#include "genotype.h"

int estimateRel(int argc, char *argv[]);
#endif
