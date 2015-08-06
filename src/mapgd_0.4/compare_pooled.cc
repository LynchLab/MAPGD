#include "compare_pooled.h"

//calculates fst between two pooled populations.
int compare_pooled(int argc, char *argv[])
{
	std::string infile="";
	std::string outfile="";

	env_t env;

	int s=0;
	real_t a=0.0;
	real_t EMLMIN=0.0001;

	std::vector <size_t> pop;

	/* Here we set all the environment options, */

	env.setname("mapgd cp");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman and Micheal Lynch");

	env.setdescription("compares allele frequencies between pooled population genomic data.");

	env.optional_arg('m',"minerror", &EMLMIN, &arg_setfloat_t, "please provide an interger", "sets minimum");
	env.optional_arg('M',"MAX", &s, &arg_setint, "please provide an interger", "sets maximum");

	env.optional_arg('p',"populations", &pop, &arg_setvectorint, "please provide a list of comma seperated integers", "choose populations to compare. Populations should be specified by comma seperated\n\t\t\t\tintigers (e.g. 1,2,3) and cannot contain spaces (e.g. 1, 2, 3 is bad). A range can\n\t\t\t\tbe specified by using a hyphen (e.g. 1-10 for populations 1 through 10) and hypenated\n\t\t\t\tranges can be strung together (e.g. 1-10,15-16) so long as ranges do not overlap.");

	env.optional_arg('a',"alpha", &a, &arg_setfloat_t, "please provide a float", "only print sites where at least one population differes from the meta population mean with a p-value less than alpha.");

	env.optional_arg('i',"in", &infile, &arg_setstr, "please provide a valid inpuit file", "specifies input file (default datain.txt)");
	env.optional_arg('o',"out", &outfile, &arg_setstr, "please provide a valid name and location for output", "specifies output file (default stdout.txt) ");

	env.flag('h',"help", &env, &flag_help, "an error occured while displaying the help mesaage", "prints this message");
	env.flag('v',"version", &env, &flag_version, "an error occured while displaying the version mesaage", "prints the program version");

	if ( parsargs(argc, argv, env) ) printUsage(env);
	if ( !env.required_set() ) printUsage(env);

	map_file in, out;		//Declare the input and output files.

	/* Open the input file. */
	if (infile.size()!=0) {
		in.open(infile.c_str(), std::fstream::in);
		if (not (in.is_open() ) ) printUsage(env);
	} else {
		in.open(std::fstream::in);
	}

	/* Open the output file. */
	if (outfile.size()!=0) {
		out.open(infile.c_str(), std::fstream::out);
		if (not (out.is_open() ) ) printUsage(env);
	} else {
		out.open(std::fstream::out);
	}

	row this_row(in.get_keys() );
	key locus_k=this_row.get_key("locus:0");
	locus *this_locus=(locus *)locus_k.begin();
	get(this_row, locus_k);	
	if (locus_k.offset()==key::nokey) 	//check to see if we found the key in the input file.
		std::cerr << __FILE__ << ":" << __LINE__ << ": no key of type locus:0 could be found. Check the column names of the input file by typeing 'mapgd view -H FILENAME'. \n";

	lnmultinomial multi(4);	

	if ( pop.size()==0 ) { 
		pop.clear();
		for (size_t x=0; x<this_locus->quartet_size(); ++x) pop.push_back(x);
	};

	while ( read_row(in, this_row)!=EOF ){
		this_locus->maskall();		//maskall quartets
		this_locus->unmask(pop);		//unmask quartets in pop
		this_locus->mask_low_cov(1);		//mask low coverage
		this_locus->sort();			//

		/* Calculate the likelihoods under the reduced model assuming no variation between populations. */

		//for (size_t x=0; x<pop.size(); ++x) llhoodP[x]=estimate_pooled(line);

		/* 4) CALCULATE THE LIKELIHOOD UNDER THE ASSUMPTION OF POPULATION SUBDIVISION. */ 

		//for (size_t x=0; x<pop.size(); ++x) { };
		
		/* Likelihood ratio test statistic; asymptotically chi-square distributed with one degree of freedom. */
		/*
		for (size_t x=0; x<pop.size(); ++x){
			if (!line.get_quartet(x).masked){
				llhoodSS=0;
				llhoodPS=0;
				for (size_t y=0; y<pop.size(); ++y){
					llhoodSS+=llhoodS[y];
					if (x!=y) llhoodPS+=llhoodS[y];
				else llhoodPS+=llhoodP[y];
				};
				llstat[x] = fabs(2.0 * (llhoodSS - llhoodPS) );
				if (llstat[x]>maxll) maxll=llstat[x];
			}
		};
		if (std::max(maxll, real_t(0) )>=a){
			for (size_t x=0; x<pop.size(); ++x){
		};*/
	}
	in.close();
	out.close();
	exit(0);
};
