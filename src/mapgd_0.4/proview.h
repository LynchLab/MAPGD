#ifndef _PROVIEW_H_
#define _PROVIEW_H_

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
#include "sample_name.h"
#include "bcf2pro-file.h"
#include "region.h"

#define SQRT2	1.41421356237

/// Struct to pass around command line arguments. Will either be depricated, or will be integreated into the interface.h*/
struct Args{
	float_t pvalue;
	int min;
	bool pro;
	bool binary;
};

/** \ingroup COMMANDS
 *  @{
 */
/// executes the proview command.
int proview(int argc, char *argv[]);
/** @}*/
#endif
