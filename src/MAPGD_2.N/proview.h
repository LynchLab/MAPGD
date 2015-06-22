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

#define SQRT2	1.41421356237

struct Args{
	float_t pvalue;
	int min;
	bool pro;
	bool binary;
};

int proview(int argc, char *argv[]);
void readheader(int *dic, Args, std::istream*, profile &);
int readline(int *dic, Args, std::istream*, profile &);
void scanCol(const std::string &, const int *, quartet_t &, const char &);

#endif
