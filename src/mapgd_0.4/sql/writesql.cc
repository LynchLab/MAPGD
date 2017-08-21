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

	if (file.indexed() && !file.concatenated()){
		std::cerr << "index?..." << std::endl;
		if (db_check_schema(db, &index)) 
		{
			std::stringstream stream;
			db_open_table(db, &index, &stream);
			while (db_get(&stream, &index) );
		}
		else 
		{
			fprintf(stderr, gettext("mapgd:%s:%d: Cannot write data because there is no SCAFFOLD table in %s. Try writing an .idx file first.\n"), __FILE__, __LINE__, db_name.c_str());
			sqlite3_close(db);
			return 0;	
		}
	}

	while(file.is_open())
	{
		std::cerr << "file is open..." << std::endl;
		db_begin(db);
		if (line->get_table_name()==File_index::table_name) 
		{
			std::cerr << "line?..." << std::endl;
			if (!db_check_schema(db, &index)) 
			{
				db_make_table(db, line);
				std::cerr << "schema..." << std::endl;
				file.read(&index);
				while (file.table_is_open() )
				{
					std::cerr << line << std::endl;
					db_insert(db, &index);
					file.read(&index);
				}
			} else {
				std::cerr << "not schema..." << std::endl;
				file.read(&index);
				while (file.table_is_open() ) file.read(&index);
			} 
		} else {
			db_make_table(db, line);
			std::cerr << "table name..." << std::endl;
			if (file.indexed() ){
				Indexed_data *indexed_line=dynamic_cast <Indexed_data *> (line);
				file.read(index, indexed_line);
				while (file.table_is_open() ){
					std::cerr << line << std::endl;
					db_insert(db, line);
					file.read(index, indexed_line);
				}
			} else {
				file.read(line);
				while (file.table_is_open() ){
					std::cerr << line << std::endl;
					db_insert(db, line);
					file.read(line);
				}
			}
		}
		db_end(db);
		delete line;
		line=file.read_header();
	}
	sqlite3_close(db);
	return 0;					//Since everything worked, return 0!.
}
