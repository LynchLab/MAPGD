#ifndef _DBAPI_
#define _DBAPI_

#include <cstdio>
#include "sqlite3.h"
#include "../typedef.h"
#include "../data_types/data.h"

int call_back(void *, int, char **, char **);

void db_begin(sqlite3 *);
void db_end(sqlite3 *);

//For writing data
void db_make_table(sqlite3 *, const Data *);
void db_insert(sqlite3 *, const Data *);

//For reading data
void db_open_table(sqlite3 *, Data *, std::stringstream *);
int db_get(std::stringstream *, Data *);


#endif

