/******************* proview.cpp *******************
 * This is inspired by sam2pro, some of the libraries
 * by Bernhard Haubold, haubold@evolbio.mpg.de
 * 			and
 * by Matthew Ackerman, matthew.s.ackerman@gmail.com
 * Description: Convert sam output to profiles.
 * Date: Wed Jul 17 12:00:00 2015
 **************************************************/

//TODO: change merging behavior to merge samples.

#include "fastview.h"

//class Bcf2pro : public Locus{};
//class Bcf2pro_file : public Indexed_file <Bcf2pro>{};

int fastview(int argc, char *argv[])
{
	/* commands that can be set from the command line */

	/* default values for arguments */

	bool out_binary=false;
	bool in_pro=false;

	std::string namefile="";
	std::string outfile="";
	std::string headerfile="";

	Environment env;
	env.set_name("mapgd fastview");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("prints data in the '.pro' file quartet format");
	env.required_arg('H',"header",	headerfile,	"You must specify an index file (-H)", "sets the index file (required to use mpileup)");
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message", "prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occurred while displaying the version message", "prints the program version");
	env.flag(	'b',"binary",	&out_binary, 	&flag_set, 	"an error occurred", "output in a binary format");
	env.flag(	'p',"pro",	&in_pro, 	&flag_set, 	"an error occurred", "input is in pro format");

	if (parsargs(argc, argv, env)!=0) exit(0);

	Indexed_file <Locus> out_file;			//the output profile
	Indexed_file <Locus> in_file;			//the output profile
	Locus in_locus;
	Locus out_locus;

	File_index index;

	in_file.open(std::ios::in);
	in_locus=in_file.read_header();

	if (outfile.size()==0) out_file.open( out_binary ? std::ios::out | std::ios::binary : std::ios::out );
	else out_file.open(outfile.c_str(), out_binary ? std::ios::out | std::ios::binary : std::ios::out );

	out_file.set_index(in_file.get_index() );
	out_file.write_header(out_locus);

	while (in_file.read(in_locus).table_is_open() ){
		out_locus=in_locus;
		out_locus.unmaskall();
		out_file.write(out_locus);
	};
	in_file.close();
	out_file.close();
	env.close();
	return 0;
}


