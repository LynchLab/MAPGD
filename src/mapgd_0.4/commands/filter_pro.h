/* A command to filter data by fields set in the map file. */

#ifndef FILTER_PRO_H_
#define FILTER_PRO_H_

#include <string>
#include "map_file.h"
#include "interface.h"
#include "data_types/allele.h"
#include "typedef.h"

#ifndef NOGSL
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>

#endif

/** 
  * \ingroup COMMANDS
  * @{
  */
/// Filters the output of the pro command.
int filter_pro(int argc, char *argv[]);
/** @} */

#endif
