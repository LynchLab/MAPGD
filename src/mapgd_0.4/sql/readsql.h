#ifndef _WRITE_SQL_H_
#define _WRITE_SQL_H_

#include "../interface.h"
#include "db_api.h"
#include "../datatypes.h"
#include "../map-file.h"

//! Reads Data classes from an SQL database
/*! @ingroup Command 
 *  Reads Data classes from an SQL database.
 */
int readsql(int argc, char *argv[]);

#endif
