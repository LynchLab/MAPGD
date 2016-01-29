#ifndef _DBAPI_
#define _DBAPI_
#include "sqlite3.h"

static int callback(void *NotUsed, int argc, char **argv, char **azColName);
void * db_begin(sqlite3 *db);
void * db_end(sqlite3 *db);
void * db_insert(sqlite3 *db, const row &this_row);
#endif
