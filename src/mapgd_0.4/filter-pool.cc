/* 

command filter:

*/

#include "filter-pool.h"

int filter_pool(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string in_file="", out_file="";

	bool binary=false;

	float_t f_min=-2., f_max=2., min_efc=0.0, error_max=0.1, max_polyll=FLT_MAX, min_polyll=-1, max_fixedll=FLT_MAX, min_fixedll=-1, min_freq=0, max_freq=1;
	int max_coverage=CNT_MAX, min_coverage=4;

	Environment env;
	env.set_name("polgd filter");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman");
	env.set_description("Filter sites in '.pol' files based on criteria.");
	env.optional_arg('E',"max-error", error_max, "please provide a real number.", "maximum error of a site (default 0.1).");
	
	env.optional_arg('i',"input", 	in_file, "please provide a filename.", "the name of the input file (default stdin).");
	env.optional_arg('o',"output", 	out_file, "please provide a filename.", "the name of an output file (default stdout).");

	env.optional_arg('c',"min-cov", min_coverage, "please provide an integer.", "minimum coverage for a population at a site for the site to be used (default 4).");
	env.optional_arg('C',"max-cov", max_coverage, "please provide an int.", "max coverage for a population at a site for the site to be used (default CNT_MAX).");
	env.optional_arg('p',"min-poly", min_polyll, "please provide a real number.", "minimum log likelihood of polymorphism (default 0.0).");
	env.optional_arg('P',"max-poly", max_polyll,"please provide a real number.", "maximum log likelihood of polymorphism (default none).");
	env.optional_arg('f',"min-hwe", min_fixedll, "please provide a real number.", "minimum log likelihood of fixed difference from reference (default 0.0).");
	env.optional_arg('F',"max-hwe", max_fixedll, "please provide a real number.", "maximum log likelihood of fixed difference from reference (default none).");
	env.optional_arg('q',"min-maf", min_freq, "please provide a real number.", "minimum frequency of minor allele (default 0.0).");
	env.optional_arg('Q',"max-maf", max_freq, "please provide a real number.", "maximum frequency of minor allele (default 0.5).");
	env.flag(	'b',"binary", 	&binary,	&flag_set, 	"please provide an int.", "output in binary mode (fast).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	float_t polyll, hwell;
	Pooled_data s;
	Indexed_file <Pooled_data> pol_in, pol_out;
	
	if (in_file.size()==0)	pol_in.open(std::fstream::in);
	else pol_in.open(in_file.c_str(), std::fstream::in);

	if (out_file.size()==0)	pol_out.open(std::fstream::out);
	else pol_out.open(out_file.c_str(), std::fstream::out);

	s=pol_in.read_header();

	pol_out.set_index(pol_in.get_index() );
	std::cerr << s.header();
	pol_out.write_header(s);

	pol_in.read(s);

	while( pol_in.table_is_open() ){
		if (s.coverage==0){
			if (min_coverage==0) pol_out.write(s);
		} else {
			if ( *(std::max_element( s.polyll.begin(),  s.polyll.end() ) ) >= min_polyll &&  *(std::max_element( s.polyll.begin(),  s.polyll.end())) <= max_polyll){ 
			if ( *(std::max_element(s.fixedll.begin(), s.fixedll.end() ) ) >= min_fixedll && *(std::max_element(s.fixedll.begin(), s.fixedll.end())) <= max_fixedll){
			if (s.error <= error_max){
			if ( *(std::min_element(s.p.begin(), s.p.end() )) >= min_freq && *(std::max_element(s.p.begin(), s.p.end() )) <= max_freq){
			if (s.coverage >= min_coverage && s.coverage <= max_coverage){
			pol_out.write(s);
		}}}}}}
		pol_in.read(s);
	}
	pol_out.close();
	return 0;					//Since everything worked, return 0!.
}
