/* 

command filter:

*/

#include "filter.h"

int filter(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string in_file="", out_file="";

	bool binary=false;

	float_t f_min=-2., f_max=2., min_efc=0.0, error_max=0.1, max_polyll=FLT_MAX, min_polyll=-1, max_hwell=FLT_MAX, min_hwell=-1, min_freq=0, max_freq=0.5, min_gof=-2.;
	count_t max_coverage=CNT_MAX, min_coverage=4, max_pitch=1;

	env_t env;
	env.setname("mapgd filter");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman");
	env.setdescription("Filter sites in '.map' files based on criteria.");
	env.optional_arg('m',"maxerror", &error_max, 	&arg_setfloat_t, "please provide a float.", "maximum error of a site (default 0.1).");
	
	env.optional_arg('i',"input", 	&in_file, 	&arg_setstr, 	"please provide an int.", "minimum coverage for an individual at a site for an individual to be used (defualt 4).");
	env.optional_arg('o',"output", 	&out_file, 	&arg_setstr, 	"please provide an int.", "the name of an output file (default std::cout).");

	env.optional_arg('c',"mincoverage", &min_coverage, 	&arg_setint, 	"please provide an int.", "minimum coverage for an individual at a site for an individual to be used (defualt 4).");
	env.optional_arg('C',"mincoverage", &max_coverage, 	&arg_setint,	"please provide an int.", "minimum coverage for an individual at a site for an individual to be used (defualt 4).");
	env.optional_arg('p',"minpoly", &min_polyll, 	&arg_setfloat_t, "please provide a float.", "minimum log likelihood of polymorphism (default 0.0).");
	env.optional_arg('P',"maxpoly", &max_polyll, 	&arg_setfloat_t, "please provide a float.", "maximum log likelihood of polymorphism (default none).");
	env.optional_arg('f',"minhwe", 	&min_hwell,	&arg_setfloat_t, "please provide a float.", "minimum log likelihood of Hardy-Weinberg disequilibrium (default 0.0).");
	env.optional_arg('F',"maxhwe", 	&max_hwell,	&arg_setfloat_t, "please provide a float.", "maximum log likelihood of Hardy-Weinberg disequilibrium (defalt none).");
	env.optional_arg('q',"minminor", &min_freq, 	&arg_setfloat_t, "please provide a float.", "minimum frequency of minor allele (default 0.0).");
	env.optional_arg('Q',"maxminor", &max_freq, 	&arg_setfloat_t, "please provide a float.", "maximum frequency of minor allele (default 0.5).");
	env.optional_arg('g',"goodfit", &min_gof,	&arg_setfloat_t, "please provide a float.", "maximum -log (p-value) that a site is bad---not corrected for multiple test (defaults 2.0).");
	env.optional_arg('N',"number", 	&max_pitch,	&arg_setint, 	"please provide an int.", "maximum number of individuals rejected by 'mapgd allele' (default 0.0).");
	env.flag(	'b',"binary", 	&binary,	&flag_set, 	"please provide an int.", "output in binary mode (fast).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occured while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.

	float_t polyll, hwell;
	allele_stat s;
	Indexed_file <allele_stat> map_in, map_out;

	map_in.open(std::fstream::in);
	map_out.open(std::fstream::out);

	s=map_in.read_header();

	map_out.set_index(map_in.get_index() );
	map_out.write_header(s);

	map_in.read(s);

	while( !map_in.eof() ){
		hwell=(s.ll-s.hwell)*2, polyll=(s.ll-s.monoll)*2;
		float_t m=s.mm+s.Mm/2.;
		if (polyll < max_polyll && hwell < max_hwell && s.error < error_max && s.f > f_min && s.f < f_max && m > min_freq && m < max_freq && polyll > min_polyll && hwell > min_hwell && s.gof > min_gof && s.efc > min_efc && s.coverage > min_coverage && s.coverage < max_coverage && s.excluded < max_pitch) map_out.write(s);
		map_in.read(s);
	}
	map_out.close();
	return 0;					//Since everything worked, return 0!.
}
