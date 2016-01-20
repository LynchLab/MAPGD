/* 

command filter:

*/

#include "filter.h"

int filter(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	bool verbose=false, binary=false;

	float_t f_min=-1., f_max=1., min_efc=0.0, error_max=0.1, min_polyll=0, min_hwell=0, min_freq=0, max_freq=0.5, min_gof=-2.;
	count_t max_coverage=MAX_COUNT, min_coverage=4, max_pitch=1;

	env_t env;
	env.setname("mapgd filter");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman");
	env.setdescription("Filter sites in map files based on criteria.");
	env.optional_arg('m',"maxerror", &error_max, 	&arg_setfloat_t, "please provide a float.", "maximum error of a site (default 0.1).");

	env.optional_arg('c',"mincoverage", &min_coverage, 	&arg_setfloat_t, 	"please provide an int.", "minimum coverage for an individual at a site for an individual to be used (defualt 4).");
	env.optional_arg('C',"mincoverage", &max_coverage, 	&arg_setfloat_t,	"please provide an int.", "minimum coverage for an individual at a site for an individual to be used (defualt 4).");
	env.optional_arg('P',"polyll", 	&min_polyll, 	&arg_setfloat_t, "please provide a float.", "cut-off value for printing polymorphic sites (default 0.0).");
	env.optional_arg('H',"hwell", 	&min_hwell,	&arg_setfloat_t, "please provide a float.", "cut-off value for printing polymorphic sites (default 0.0).");
	env.optional_arg('f',"minfreq", &min_freq, 	&arg_setfloat_t, "please provide a float.", "cut-off value for printing polymorphic sites (default 0.0).");
	env.optional_arg('F',"maxfreq", &max_freq, 	&arg_setfloat_t, "please provide a float.", "cut-off value for printing polymorphic sites (default 0.0).");
	env.optional_arg('g',"goodfit", &min_gof,	&arg_setfloat_t, "please provide a float.", "cut-off value for the goodness of fit statistic (defaults 2.0).");
	env.optional_arg('N',"number", 	&max_pitch,	&arg_setint, 	"please provide an int.", "cut-off value for number of bad individuals accepted at a site (default "+str(max_pitch)+").");
	env.optional_arg( , , 		&max_pitch,	&arg_setint, 	"please provide an int.", "cut-off value for number of bad individuals needed before a site is removed entirely (default 96).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message");
	env.flag(	'b',"binary", 	&binary,	&flag_set, 	"please provide an int.", "input is binary");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occured while displaying the version message.", "prints the program version");
//	env.flag(	'V', "verbose", &verbose,	&flag_set, 	"an error occured while enabeling verbose excecution.", "prints more information while the command is running.");
//	env.flag(	'q', "quiet", 	&quite,		&flag_set, 	"an error occured while enabeling quite execution.", "prints less information while the command is running.");

	if ( parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.


	alelel_stat s;
	char major, minor, ref;
	std::string id1, id2;
	in >> id1 >> id2 >> ref >> major >> minor >> s;
	float_t hwell=s.hwell-s.monoll, polyll=s.polyll-s.monoll;
	if (s.error < error_max && s.f > f_min && s.f < f_max && m > m_min && m < m_max && polyll > min_polyll && hwell > min_hwell && s.gof > min_gof && s.efc > min_efc && s.coverage > min_coverage && s.coverage < max_coverage && s.excluded < max_pitch)
	out << ids << ids << major << minor << s;
	
	return 0;					//Since everything worked, return 0!.
}
