/* 

command filter:

*/

#include "readsql.h"

int readsql(int argc, char *argv[])
{

	std::string db_name="";
	std::string query="";
	std::string table="";
	/* All the variables that can be set from the command line */

	Environment env;
	env.set_name("mapgd read");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("reads mapgd output from an sql database.");
	env.required_arg('d',"database", db_name,"please the name of a database", "the name of the database");
	env.optional_arg('q',"query", 	 query,"please provide an str", "runs a query on the database.");
	env.optional_arg('t',"table", 	 table,"please provide an str", "prints the table.");
	env.flag(	'v', "version",  &env, 		&flag_version, 	"an error occurred while displaying the version message", "prints the program version");
	env.flag(	'h', "help", 	 &env, 		&flag_help, 	"an error occurred while displaying the version message", "prints the program version");

	if (parsargs(argc, argv, env)!=0) print_usage(env); //Gets all the command line options, and prints usage on failure.

	sqlite3 *db;
	int rc;

	rc = sqlite3_open(db_name.c_str(), &db);
	if( rc ){
		fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
		sqlite3_close(db);
		return(1);
	} else {
		fprintf(stderr, "Opened database successfully\n");
	}

 	Base_file file;
 	File_index index;
 	file.open(std::ios::out);
	std::vector <std::string> columns={};

	if (table!="")
	{
		fprintf(stderr, "Trying to open data of type %s\n", table.c_str());
		Data *data=Data::new_from_str(table, columns);
		fprintf(stderr, "done.\n");
		fprintf(stderr, "New data of type %s\n", data->get_table_name().c_str() );
		fprintf(stderr, "New data of type %s\n", data->table_name.c_str() );
		fprintf(stderr, "Print once? %d\n", data->get_print_once() );
		columns=db_get_constructor(db, data);
		delete data;
		data=Data::new_from_str(table,columns);
		std::stringstream stream;
		if (data->indexed() )
		{
			db_open_table(db, &index, &stream);
			while (db_get(&stream, &index) );
			stream.clear();

			Indexed_data *indexed_data=dynamic_cast <Indexed_data *> (data);

			if (query!="")
				db_open_table_w_query(db, indexed_data, &stream, query);
			else
				db_open_table(db, indexed_data, &stream);

			file.write_header(index, indexed_data);
			while (db_get(&stream, indexed_data) ) 
				if (!indexed_data->get_print_once() )
					file.write(index, indexed_data);
			file.write(index, indexed_data);
			file.close();
		}
		else 
		{
			if (query!="")
				db_open_table_w_query(db, data, &stream, query);
			else
				db_open_table(db, data, &stream);
			file.write_header(data);
			while (db_get(&stream, data) ) 
				if (!data->get_print_once() )
					file.write(data);
			file.write(data);
			file.close();
		}
	}
	return 0;					//Since everything worked, return 0!.
}
