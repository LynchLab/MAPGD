/***** proview.h **********************************
 * Description: Convert mpileup output to pro-files.
 */

#ifndef PROVIEW_H_
#define PROVIEW_H_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "interface.h"
#include "stream-tools.h"
#include "pro-file.h"

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
