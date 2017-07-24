/* 

command filter:

*/

#include "test_keys.h"

int test_keys(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	bool all=false;
	std::string quarry="";

	Environment env;
	env.set_name("mapgd keyinfo");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Test all registered data types to ensure they are defined.");
	
        env.positional_arg('k',"key", quarry,        "No key specified", "a key to be described");

	env.flag(	'a', "avail", 	&all, 		&flag_set, 	"an error occurred while displaying the help message.", "list all known keys");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	if (quarry=="" && !all) print_usage(env);

	Flat_file <Key> keys;  	
	Key key;
	std::string filename(std::string(PATH)+"../keys/keys.txt");
	keys.open(filename.c_str(), READ);
	std::map <std::string, Key *> keymap;
	std::vector <std::string> empty, column_names;
	std::stringstream ss;

	keys.read_header();
	while (keys.table_is_open()){
		keys.read(key);
		keymap[std::string(key.name)]=new Key(key);
		if (all) std::cout << key.name << std::endl;
	}
	if (all) return 0;
	std::cout << key.header();
	std::map <std::string, Key *>::iterator ret=keymap.find(quarry);	
	if (ret==keymap.end() ) std::cerr << "Could not find key. Tell Matt to handle fuzzy matching." << std::endl;
	else std::cout << *(ret->second) << std::endl;

/*	std::vector<std::string> key_names=registry_list();	
	for ( std::vector <std::string> ::iterator it=key_names.begin(); it!=key_names.end(); it++) 
	{
		std::cerr << *it << std::endl;
		ss << Data::new_from_str(*it, empty)->header() << std::endl;
		column_names=split(ss, keys.get_delim() );
		for ( std::vector <std::string> ::iterator jt=column_names.begin(); jt!=column_names.end(); jt++) 
		{
			int start=0, end=jt->size();
			if((*jt)[start]=='@') start++;
			if((*jt)[end-1]=='\n') end--;
			std::string trimmed=jt->substr(start, end);
			std::map <std::string, Key *>::iterator ret=keymap.find(trimmed);
		//	if (ret==keymap.end() ) std::cerr << trimmed << " not found\n";
		//	else std::cerr << trimmed << " found\n";
		}
	}*/

	return 0;					//Since everything worked, return 0!.
}
