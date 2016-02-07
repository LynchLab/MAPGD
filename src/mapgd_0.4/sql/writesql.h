#ifndef _WRITE_SQL_H_
#define _WRITE_SQL_H_

#include "../interface.h"
#include "db_api.h"
#include "../datatypes.h"
#include "../map-file.h"

//! Writes Data classes to an SQL database
/*! @ingroup Command 
 *  Writes Data classes to an SQL database.
 */
int writesql(int argc, char *argv[]);

#endif
