/* 

command filter:

*/

#include "readsql.h"

int readsql(int argc, char *argv[])
{

	std::string db_name="";
	std::string query="";
	/* All the variables that can be set from the command line */

	env_t env;
	env.set_name("mapgd read");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("reads mapgd output from an sql database.");
	env.required_arg('d',"database", &db_name, 	&arg_setstr, 	"please the name of a database", "the name of the database");
	env.required_arg('q',"query", 	 &query, 	&arg_setstr, 	"please provide an str", "query");
	env.flag(	'v', "version",  &env, 		&flag_version, 	"an error occured while displaying the version message", "prints the program version");
	env.flag(	'h', "help", 	 &env, 		&flag_help, 	"an error occured while displaying the version message", "prints the program version");

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

	Indexed_file <Locus> locus_file;
	Flat_file <File_index> index_file;
	File_index index;
	std::stringstream stream;
	db_open_table(db, &index, &stream);
	while (db_get(&stream, &index) );
	stream.clear();
	Locus locus;
	db_open_table(db, &locus, &stream);

	locus_file.open(std::ios::out);
	locus_file.set_index(index);
	locus_file.write_header(locus);

	while (db_get(&stream, &locus) ) locus_file.write(locus);

	locus_file.close();
	return 0;					//Since everything worked, return 0!.
}
