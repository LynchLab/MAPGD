/* 

command filter:

*/

#include "filter.h"

int filter(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string in_file="", out_file="";

	bool binary=false;

	float_t f_min=-2., f_max=2., min_efc=0.0, error_max=0.1, max_polyll=FLT_MAX, min_polyll=-1, max_hwell=FLT_MAX, min_hwell=-1, min_freq=0, max_freq=0.5, min_gof=2., pbias=0.0;
	int max_coverage=CNT_MAX, min_coverage=4, max_pitch=1;

	Environment env;
	env.set_name("mapgd filter");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Filter sites in '.map' files based on criteria.");
	env.optional_arg('E',"max-error", error_max, "please provide a real number.", "maximum error of a site (default 0.1).");
	
	env.optional_arg('i',"input", 	in_file, "please provide a filename.", "the name of the input file (default stdin).");
	env.optional_arg('o',"output", 	out_file, "please provide a filename.", "the name of an output file (default stdout).");

	env.optional_arg('c',"min-cov", min_coverage, "please provide an integer.", "minimum coverage for a population at a site for the site to be used (default 4).");
	env.optional_arg('C',"max-cov", max_coverage, "please provide an int.", "max coverage for a population at a site for the site to be used (default CNT_MAX).");
	env.optional_arg('p',"min-poly", min_polyll, "please provide a real number.", "minimum log likelihood of polymorphism (default 0.0).");
	env.optional_arg('P',"max-poly", max_polyll,"please provide a real number.", "maximum log likelihood of polymorphism (default none).");
	env.optional_arg('f',"min-hwe", min_hwell, "please provide a real number.", "minimum log likelihood of Hardy-Weinberg disequilibrium (default 0.0).");
	env.optional_arg('F',"max-hwe", max_hwell, "please provide a real number.", "maximum log likelihood of Hardy-Weinberg disequilibrium (default none).");
	env.optional_arg('X',"pbias",   pbias, "please provide a real number.", "maximum likelihood of reference bias (default none).");
	env.optional_arg('q',"min-maf", min_freq, "please provide a real number.", "minimum frequency of minor allele (default 0.0).");
	env.optional_arg('Q',"max-maf", max_freq, "please provide a real number.", "maximum frequency of minor allele (default 0.5).");
	env.optional_arg('g',"min-fit", min_gof, "please provide a real number.", "maximum -log (p-value) that a site is bad---not corrected for multiple test (defaults 2.0).");
	env.optional_arg('N',"number", 	max_pitch, "please provide an integer.", "maximum number of individuals rejected by 'mapgd allele' (default 0.0).");
	env.flag(	'b',"binary", 	&binary,	&flag_set, 	"please provide an int.", "output in binary mode (fast).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");
	
	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	float_t polyll, hwell;
	Allele s;
	Indexed_file <Allele> map_in, map_out;
	
	if (in_file.size()==0)	map_in.open(std::fstream::in);
	else map_in.open(in_file.c_str(), std::fstream::in);

	if (out_file.size()==0)	map_out.open(std::fstream::out);
	else map_out.open(out_file.c_str(), std::fstream::out);

	s=map_in.read_header();

//	std::cerr << "DEL:" << s.delim << ".\n";

	map_out.set_index(map_in.get_index() );
	map_out.write_header(s);

	map_in.read(s);

	while( map_in.table_is_open() ){
		hwell=(s.ll-s.hwell)*2, polyll=(s.ll-s.monoll)*2;
		float_t m=s.mm+s.Mm/2.;
		if (s.coverage==0){
			if (min_coverage==0) map_out.write(s);
		} else {
			if (polyll >= min_polyll){
			if (hwell >= min_hwell && hwell <= max_hwell){
			if (s.error <= error_max){
			if (s.f >= f_min && s.f <= f_max){
			if (m >= min_freq && m <= max_freq){
			if (s.gof >= -min_gof){
			if (s.efc >= min_efc){
			if (s.coverage >= min_coverage && s.coverage <= max_coverage){
			if (s.excluded <= max_pitch)
			if (s.pbias >= pbias)
				map_out.write(s);
		}}}}}}}}}
		map_in.read(s);
	}
	map_out.close();
	return 0;					//Since everything worked, return 0!.
}
