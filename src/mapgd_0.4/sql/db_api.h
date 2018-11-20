#ifndef _DBAPI_
#define _DBAPI_

#include <cstdio>
#include "sqlite3.h"
#include "../typedef.h"
#include "../data_types/data.h"
#include "stream_tools.h"

int call_back(void *, int, char **, char **);

void db_begin(sqlite3 *);
void db_end(sqlite3 *);

//!For writing data
void db_make_table(sqlite3 *, const Data *);
void db_insert(sqlite3 *, const Data *);

//!For reading data
void db_open_table(sqlite3 *, Data *, std::stringstream *);
void db_open_table_w_query(sqlite3 *, Data *, std::stringstream *, const std::string &);
int db_get(std::stringstream *, Data *);

std::vector <std::string> db_get_constructor(sqlite3 *, Data *);

int db_check_schema(sqlite3 *, Data *);

#endif

