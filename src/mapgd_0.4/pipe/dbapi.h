static int callback(void *NotUsed, int argc, char **argv, char **azColName){
   int i;
   for(i=0; i<argc; i++){
      printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
   }
   printf("\n");
   return 0;
}


void * db_begin(sqlite3 *db)
{
	sqlite3_exec(db, "BEGIN TRANSACTION;", callback, 0, &error_message);
};

void * db_end(sqlite3 *db)
{
        sqlite3_exec(db, "END TRANSACTION;", callback, 0, &error_message);
};

void * db_insert(sqlite3 *db, const row &this_row)
{
	char make_table1[255]={0};
	char make_table2[255]={0};
	char add_data1[255]={0};
	char add_data2[255]={0};
	snprintf (make_table1, 255, "CREATE TABLE IF NOT EXISTS %s (\n", this_row.table() );
	snprintf (add_data1, 255, "INSERT INTO %s VALUES (", this_row.table() );
	for ( auto &this_key : this_row.keys_){
		if (&this_key!=&this_row.keys_.back() ) {
			snprintf (add_data2, 255, "%s %s,", add_data1, this_key.to_string( fetch(this_row, this_key) ) );
			snprintf (make_table2, 255, "%s\t%s\t%s\n,", make_table1, this_key.name(), this_key.sql_type() );
		}
		else{
			snprintf(add_data2, 255, "%s %s", add_data1, this_key.to_string( fetch(this_row, this_key) ) );
			snprintf (make_table2, 255, "%s\t%s\t%s\n", make_table1, this_key.name(), this_key.sql_type() );
		}
		memcpy(make_table1, make_table2, 255);
		memcpy(add_data1, add_data2, 255);
	};
	snprintf (add_data2, 255, "%s );", add_data1 );
	snprintf (make_table2, 255, "%s );", make_table1 );
	sqlite3_exec(db, make_table2, callback, 0, &error_message);
	sqlite3_exec(db, add_data2, callback, 0, &error_message);
}

