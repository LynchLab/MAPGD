#ifndef _READ_SQL_H_
#define _READ_SQL_H_

#include <libintl.h>
#include "interface.h"
#include "db_api.h"
#include "datatypes.h"
#include "map_file.h"

//! Reads Data classes from an SQL database
/*! @ingroup Command 
 *  Reads Data classes from an SQL database.
 */
int readsql(int argc, char *argv[]);

#endif
