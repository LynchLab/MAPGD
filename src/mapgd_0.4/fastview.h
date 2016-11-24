#ifndef _FASTVIEW_H_
#define _FASTVIEW_H_

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

#define SQRT2	1.41421356237

/** \ingroup COMMANDS
 *  @{
 */
/// executes the proview command.
int fastview(int argc, char *argv[]);
/** @}*/
#endif
