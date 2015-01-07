/***** sam2pro.h **********************************
 * Description: Convert sam output to profiles.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jul 21 22:46:11 2010
 **************************************************/

#ifndef PROVIEW_H_
#define PROVIEW_H_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "interface.h"
#include "streamtools.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/math/special_functions/gamma.hpp>


#define SQRT2	1.41421356237


struct Args{
	double pvalue;
	int min;
	int c;
	char delim;
	int notrim;
};

int proview(int argc, char *argv[]);
void runAnalysis(int *dic, Args, std::istream*, std::ostream*);
void scanCol(std::string column, int *dic, int *count, char consensus);

#endif
