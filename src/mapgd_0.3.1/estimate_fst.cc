/* 
*/

#include "estimate_individual.h"

int estimate_fst(int argc, char *argv[])
{
	std::vector <int> pop1={0}, pop2={1};

	/* sets up the help messages and options, see the 'interface.h' for more detials. */

	env_t env;
	env.setname("mapgd fst");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman and Takahiro Maruki");
	env.setdescription("Uses a maximum likelihood approach to estimate fst between population.");

	env.required_arg('1',"individuals", &pop1, 	&arg_setvectorint, "please provide a list of integers", "the individuals in population 1.\n\t\t\t\ta comma seperated list containing no spaces, and the format X-Y can be used to specify a range (defualt 1).");
	env.required_arg('2',"individuals", &pop2, 	&arg_setvectorint, "please provide a list of integers", "the individuals in population 2.\n\t\t\t\ta comma seperated list containing no spaces, and the format X-Y can be used to specify a range (defualt 2).");

	if ( parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.

	profile pro;		//profile is a fairly complete class that lets us read and write from pro files, 
				//which are files containing set of read 'quartets' that specify the number of 
				//A,C,G and T read at some specific location in a genome. See proFile.h for more info.

	std::ostream *out=&std::cout;
	pro.open(std::iostream::in);			
	std::vector <float_t> gofs(pro.size() );
	Locus line;
	models model;
	allele_stat mle1, mle2, mle_both;
	
/*	if(pro.read(line)!=EOF){
		line.maskall();
		line.unmask(&pop1);
		mle1=estimate (line, model, gofs, 0, 0, 0, 0);
		line.maskall();
		line.unmask(&pop2);
		mle2=estimate (line, model, gofs, 0, 0, 0, 0);
		line.maskall();
		line.unmask(&pop1);
		line.unmask(&pop2);
		mle_both=estimate (line, model, gofs, 0, 0, 0, 0);
		if (2*(mle_both.ll-mle_both.monoll)>=A){
			*out << std::fixed << std::setprecision(6) << pro.getids(line) << '\t' << (mle1.h+mle2.h)/(2*mle_both.h) << '\t'  << mle1.ll+mle2.ll-mle_both.ll << std::endl;
		}
	}*/
	pro.close();
	return 0;					//Since everything worked, return 0!.
}
