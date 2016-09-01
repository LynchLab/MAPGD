/* A command to filter data by fields set in the map file. */

#ifndef FILTER_H_
#define FILTER_H_

#include <string>
#include "map-file.h"
#include "interface.h"
#include "data_types/allele.h"
#include "typedef.h"

/** 
  * \ingroup COMMANDS
  * @{
  */
/// Filters the output of the allele command.
int filter(int argc, char *argv[]);
/** @} */

#endif
