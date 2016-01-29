/***** proview.h **********************************
 * Description: Convert mpileup output to pro-files.
 */

#ifndef PROVIEW_H_
#define PROVIEW_H_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include "interface.h"
#include "stream-tools.h"
#include "map-file.h"
#include "locus.h"

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

#endif
