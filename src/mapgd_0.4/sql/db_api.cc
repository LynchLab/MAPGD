#include "db_api.h"

static int callback(void *NotUsed, int argc, char **argv, char **azColName){
   int i;
   for(i=0; i<argc; i++){
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
   }
   printf("\n");
   return 0;
}


void db_begin(sqlite3 *db)
{
	char *error_message = 0;
	sqlite3_exec(db, "BEGIN TRANSACTION;", callback, 0, &error_message);
}

void db_end(sqlite3 *db)
{
	char *error_message = 0;
        sqlite3_exec(db, "END TRANSACTION;", callback, 0, &error_message);
}

void db_make_table(sqlite3 *db, const Data &these_data)
{
	char *error_message = 0;
	char make_table[255]={0};
	snprintf (make_table, 255, "CREATE TABLE IF NOT EXISTS %s (%s);\n", Data::table_name.c_str(), these_data.sql_header().c_str() );
	sqlite3_exec(db, make_table, callback, 0, &error_message);
}

void db_insert(sqlite3 *db, const Data &these_data)
{
	char *error_message = 0;
	char add_data[255]={0};
	snprintf (add_data, 255, "INSERT INTO %s (%s) VALUES (%s);\n", Data::table_name.c_str(), these_data.sql_header().c_str(), these_data.sql_values().c_str() );
	sqlite3_exec(db, add_data, callback, 0, &error_message);
}
