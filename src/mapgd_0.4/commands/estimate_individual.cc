#include "estimate_individual.h"

#define BUFFER_SIZE 500			//Sets the size of the row buffer.
#define PRAGMA				//ENABLES openMP pragmas.


/// Estimates a number of summary statistics from short read sequences.
/**
 */
void estimate (locus &site, row &out, models &model, const count_t &MIN, const float_t &EMLMIN, const float_t &MINGOF, const count_t &MAXPITCH){

	allele mle, temp;					//allele is a basic structure that containes all the summary statistics for
								//an allele. It gets passed around a lot, and a may turn it into a class that has
								//some basic read and write methods.
	site.mask_low_cov(MIN);
	size_t texc=site.maskedcount(), rexc;
	std::vector <real_t> gofs(site.quartet_size(),0);
	rexc=texc;

	if (init_params(site, mle, EMLMIN) ){			//If >90% of reads agree, then assume a homozygote,
								//otherwise, assume heterozygote.
	if (site.get_count(1)!=0){
		rexc=maximize_grid(site, mle, model, gofs, -MINGOF, MAXPITCH+texc);		//trim bad clones and re-fit the model.
//		rexc=maximize_newton(site, mle, model, gofs, MIN, -MINGOF, MAXPITCH+texc);	//trim bad clones and re-fit the model.
	}
	else
		rexc=maximize_analytical(site, mle, model, gofs, -MINGOF, MAXPITCH+texc);	//trim bad clones and re-fit the model.
	}

	size_t excluded=rexc-texc;

	// CALCULATE THE LIKELIHOODS 

	real_t maxll=model.loglikelihood(site, mle);		//Sets the site.ll to the log likelihood of the best fit (ll). 

	if (get_freq(mle)<0.5){					//Check to see if the major and minor alleles are reversed.
		std::swap(mle.base[0], mle.base[1]);
		std::swap(mle.frequency[0], mle.frequency[2]);
		site.swap(0, 1);
		
	}
	else if (get_freq(mle)==0.5){
		if (rand() % 2){				//If the major and minor allele frequencies are identical, 
			std::swap(mle.base[0], mle.base[1]);	//flip a coin to determine the major and minor allele.
			std::swap(mle.frequency[0], mle.frequency[2]);
			site.swap(0, 1);
		};
	};

	temp=mle; temp.frequency[0]=1.0; 
	temp.frequency[1]=0.; temp.frequency[2]=0.;		//Copies site to mono, then sets mono to a monomophic site 
								//(i.e. sets the genotypic frequencies Mm and mm to 0.
	temp.error=(real_t)(site.get_coverage()-site.get_count(0) )/real_t(site.get_coverage() );//Sets the error rate of mono to the null error rate.
	real_t monoll=model.loglikelihood(site, temp);		//Sets the site.ll to the log likelihood of the best fit (ll). 

	if (monoll>maxll){
		maxll=monoll;
		mle.error=temp.error;
	};

	temp=mle; 

	temp.frequency[0]=pow(get_freq(mle), 2);				//Similar set up to mono, but now assuming 
	temp.frequency[1]=2.*get_freq(mle)*(1.-get_freq(mle) ); 		//Hardy-Weinberg equilibrium.
	temp.frequency[2]=pow(1.-get_freq(mle) , 2);				//?
	real_t hwell=model.loglikelihood(site, temp);			//?
	real_t this_efc=efc(site);

	real_t m1=maxll-monoll, m2=maxll-hwell;  
	place(out, polyll_k,	&m1);
	place(out, hwell_k, 	&m2);
	place(out, ncut_k, 	&excluded);
	place(out, efchrom_k, 	&this_efc );
	place(out, maxll_k, 	&maxll);
	place(out, allele_k, 	&mle);
}


/// The estimate_individual command 
int estimate_individual(int argc, char *argv[])
{
	/* Default values for the variables that can be set from the command line. I may move some of these to a
	 * configuration file of some kind, but linux programs are typically desgned to be run as simple binary files
	 * that can reside anywhere on a computer. As such, I can't really depend on being in the same directory as a
	 * configuration file most of the time. I may make a .rc file or something that will overide these values if it 
	 * is present. 
         */

	std::string infile="";
	std::string outfile="";

//	table in, out, env;

	bool verbose=false;
	bool quite=false;
	bool noheader=false;

	float_t EMLMIN=0.001;
	count_t MIN=0;
	float_t A=0.00;
	float_t MINGOF=2.00;
	count_t MAXPITCH=96;

	std::vector <int> ind;

	/* sets up the help messages and options, see the 'interface.h' for more detials. */

	env_t env;
	env.setname("mapgd estimate");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman and Takahiro Maruki");
	env.setdescription("Uses a maximum likelihood approach to estimate population genomic statistics from a population.");
	env.optional_arg('I', "individuals", &ind, 	&arg_setvectorint, "please provide a list of integers", "the individuals to be used in estimates.\n\t\t\t\ta comma seperated list containing no spaces, and the format X-Y can be used to specify a range (defualt ALL).");
	env.optional_arg('m', "minerror", &EMLMIN, 	&arg_setfloat_t, "please provide a float.", "prior estimate of the error rate (defualt 0.001).");
	env.optional_arg('M', "mincoverage", &MIN, 	&arg_setint, 	"please provide an int.", "minimum coverage of sites to be estimated (defualt 4).");
	env.optional_arg('a', "alpha", 	&A, 		&arg_setfloat_t, "please provide a float.", "cut-off value for printing polymorphic sites (default 0.0).");
	env.optional_arg('g', "goodfit", &MINGOF,	&arg_setfloat_t, "please provide a float.", "cut-off value for the goodness of fit statistic (defaults 2.0).");
	env.optional_arg('N', "number", 	&MAXPITCH,	&arg_setint, 	"please provide an int.", "maximum number of clones to be trimmed (default 96).");
	env.flag(	 'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message");
	env.flag(	 'v', "version", &env, 		&flag_version, 	"an error occured while displaying the version message.", "prints the program version");
	env.flag(	 'V', "verbose", &verbose,	&flag_set, 	"an error occured while enabeling verbose excecution.", "prints more information while the command is running.");
	env.flag(	 'q', "quite", 	&quite,		&flag_set, 	"an error occured while enabeling quite execution.", "prints less information while the command is running.");

	if ( parsargs(argc, argv, env) ) printUsage(env); 	//Gets all the command line options, and prints usage on failure.

	models model;						//model is needed to do our likelihood calculations.

	std::list <key> required={allele_k};			// makes a key to read/write a "ALLELE_STATS".


	row  in_row[BUFFER_SIZE]=row( in.get_keys("PROFILE") );
	row out_row[BUFFER_SIZE]=row( in.get_keys("GENOME").Union(required).ToList() );

	while (true){			
		size_t c=0, readed=0;
		bool estimate_me=1;
		#ifdef PRAGMA
		#pragma omp parallel private(c, model, estimate_me) 
		#endif
		{
			locus* this_locus;
			allele* this_allele;
			key locus_k=in_row[0].get_key("LOCUS");
			#ifdef PRAGMA
			#pragma omp for
			#endif
			for (size_t x=0; x<BUFFER_SIZE; ++x){
				#ifdef PRAGMA
				#pragma omp critical
				#endif
				{
					c=readed;				//Turn on the ability to read data from all clones in 
					if(read_row(in, in_row[c])!=EOF){
						readed++;			//reads the next line of the std::in. Retuerns 0 if no next line.
						estimate_me=1;
					}
					else estimate_me=0;
				}
				if(estimate_me) {
					this_locus=(locus *)fetch(in_row[c], locus_key);			//get a locus pointer with locus_key.
					estimate(this_locus, out_row[c], model, MIN, EMLMIN, MINGOF, MAXPITCH);	//estiamte me.
				}
			}

		}
                for (uint32_t x=0; x<readed; ++x){
			if (2*(bestll-monoll)>=A) write_row(out, out_row[x]);	// write the output. 
		}
		if (readed!=BUFFER_SIZE){break;}
	}
	in.close();
	out.close();
	return 0;					//Since everything worked, return 0!.
}
