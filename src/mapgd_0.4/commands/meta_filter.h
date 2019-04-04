/* A command to filter data by fields set in the map file. */

#ifndef META_FILTER_H_
#define META_FILTER_H_

#include <string>
#include "map_file.h"
#include "interface.h"
#include "data_types/allele.h"
#include "typedef.h"

/** 
  * \ingroup COMMANDS
  * @{
  */
/// Filters the output of the allele command.
int meta_filter(int argc, char *argv[]);
/** @} */

#endif
