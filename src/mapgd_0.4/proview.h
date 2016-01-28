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
#include "file-index.h"
#include "map-file.h"

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#define SQRT2	1.41421356237

/*! \breif struct to pass around command line arguments. Will either be depricated, or will be integreated into the interface.h*/
struct Args{
	float_t pvalue;
	int min;
	bool pro;
	bool binary;
};

/*! \breif executes the proview command.*/
int proview(int argc, char *argv[]);

/*! \breif Depricated?*/
void readheader(int *dic, Args, std::istream*, profile &);

/*! \breif Depricated?*/
int readline(int *dic, Args, std::istream*, profile &);

/*! \breif Bernhard's scanCol function. Nice and fast.*/
void scanCol(const std::string &, const int *, quartet_t &, const char &);

#endif
