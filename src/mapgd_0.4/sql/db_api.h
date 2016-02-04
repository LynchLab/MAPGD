#ifndef _DBAPI_
#define _DBAPI_

#include <cstdio>
#include "sqlite3.h"
#include "../data_types/data.h"

int call_back(void *, int, char **, char **);
void db_begin(sqlite3 *);
void db_end(sqlite3 *);
void db_make_table(sqlite3 *, const Data *);
void db_insert(sqlite3 *, const Data *);
#endif
