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
#include "proFile.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
//#include <boost/math/special_functions/gamma.hpp>


#define SQRT2	1.41421356237


struct Args{
	float_t pvalue;
	int min;
	int c;
	char qdel, cdel;
	bool notrim;
	bool noheader;
	bool pro;
};

int proview(int argc, char *argv[]);
void runAnalysis(int *dic, Args, std::istream*, profile &);
void scanCol(std::string column, int *dic, count_t *count, char consensus);

#endif
