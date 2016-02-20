/* 

command filter:

*/

#include "writesql.h"

int writesql(int argc, char *argv[])
{

	std::string db_name="";
	/* All the variables that can be set from the command line */

	Environment env;
	env.set_name("mapgd writesql");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("writes mapgd output to an sql database.");
	env.required_arg('d',"database", db_name,"please provide an str.", "the name of destination database");
	env.flag(	'v', "version",  &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");
	env.flag(	'h', "help", 	 &env, 		&flag_help, 	"an error occurred while displaying the version message.", "prints the program version");

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
	file.open(std::ios::in);
	Data *line=file.read_header();

	while(file.is_open() ){
		db_begin(db);
		db_make_table(db, line);
		if (line->get_table_name()==File_index::table_name) {
			file.read(&index);
			while (file.table_is_open() ){
				db_insert(db, &index);
				file.read(&index);
			}
		} else {
			if (file.indexed() ){
				Indexed_data *indexed_line=dynamic_cast <Indexed_data *> (line);
				file.read(index, indexed_line);
				while (file.table_is_open() ){
					db_insert(db, line);
					file.read(index, indexed_line);
				}
			} else {
				file.read(line);
				while (file.table_is_open() ){
					db_insert(db, line);
					file.read(line);
				}
			}
		}
		db_end(db);
		line=file.read_header();
	}
	sqlite3_close(db);
	return 0;					//Since everything worked, return 0!.
}
