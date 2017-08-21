#include "db_api.h"
#include <iostream>
#include <sstream>

int call_back(void *stream, int argc, char **argv, char **azColName){
	std::stringstream *this_stream=(std::stringstream *)stream;
	for(int i=0; i<argc; i++){
		argv[i] ? *this_stream << argv[i] << '\t' : *this_stream << "NULL\t";
	}
	*this_stream << std::endl;
	return 0;
}

void db_begin(sqlite3 *db)
{
	char *error_message = 0;
	sqlite3_exec(db, "BEGIN TRANSACTION;\n", call_back, 0, &error_message);
#ifdef DEBUG
	std::cerr << "BEGIN TRANSACTION;\n";
#endif
}

void db_end(sqlite3 *db)
{
	char *error_message = 0;
        sqlite3_exec(db, "END TRANSACTION;\n", call_back, 0, &error_message);
#ifdef DEBUG
	std::cerr << "END TRANSACTION;\n";
#endif
}

void db_make_table(sqlite3 *db, const Data *these_data)
{
	char *error_message = 0;
	char make_table[SQL_LINE_SIZE]={0};
	snprintf (make_table, SQL_LINE_SIZE, "CREATE TABLE IF NOT EXISTS %s %s;\n", these_data->get_table_name().c_str(), these_data->sql_header().c_str() );
	sqlite3_exec(db, make_table, call_back, 0, &error_message);
}

void db_insert(sqlite3 *db, const Data *these_data)
{
	char *error_message = 0;
	char add_data[SQL_LINE_SIZE]={0};
	snprintf (add_data, SQL_LINE_SIZE, "INSERT INTO %s %s VALUES %s;\n", these_data->get_table_name().c_str(), these_data->sql_column_names().c_str(), these_data->sql_values().c_str() );
	sqlite3_exec(db, add_data, call_back, 0, &error_message);
}


void db_open_table(sqlite3 *db, Data *these_data, std::stringstream *stream_ptr)
{
	char *error_message = 0;
	char get_data[SQL_LINE_SIZE]={0};
	sqlite3_stmt *query;
	//snprintf (get_data, SQL_LINE_SIZE, "%s;\n", these_data->sql_constructor().c_str(); );
	//sqlite3_exec(db, get_data, call_back, stream, &error_message);
	//these_data->set(stream);
	snprintf (get_data, SQL_LINE_SIZE, "SELECT * from %s;\n", these_data->get_table_name().c_str() );
	sqlite3_exec(db, get_data, call_back, stream_ptr, &error_message);
}

void 
db_close_table(sqlite3_stmt *query)
{
	sqlite3_finalize(query);
}

int 
db_get(std::stringstream *stream_ptr, Data *these_data)
{
	these_data->sql_read(*stream_ptr);
	return stream_ptr->peek()!=EOF;
}

std::vector <std::string>
db_get_constructor(sqlite3 *db, Data *these_data)
{
	char *error_message = 0;
	std::stringstream stream;
	char columns[SQL_LINE_SIZE]={0};
	snprintf (columns, SQL_LINE_SIZE, "SELECT DISTINCT SMPNUM FROM %s;", these_data->get_table_name().c_str() );
	sqlite3_exec(db, columns, call_back, &stream, &error_message);
	std::vector <std::string> str={"@ID0","ID1","REF"};
	std::string line;
	while (std::getline(stream, line) )
	{
		str.push_back(line);
	}
	return str;
}

int 
db_check_schema(sqlite3 *db, Data *these_data)
{
	char *error_message = 0;
	char check_table[SQL_LINE_SIZE]={0};
	snprintf (check_table, SQL_LINE_SIZE, "SELECT name FROM sqlite_master WHERE type='table' AND name='%s';", these_data->get_table_name().c_str() );
	std::stringstream stream;
	sqlite3_exec(db, check_table, call_back, &stream, &error_message);
	std::string str;
	stream >> str;
	if (str==these_data->get_table_name()) return 1;
	return 0;
}
