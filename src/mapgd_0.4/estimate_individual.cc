#include "estimate_individual.h"

#define BUFFER_SIZE 500			//Sets the size of the row buffer.
#define PRAGMA				//ENABLES openMP pragmas.


/// Estimates a number of summary statistics from short read sequences.
/**
 */
row estimate (locus &site, models &model, const count_t &MIN, const float_t &EMLMIN, const float_t &MINGOF, const count_t &MAXPITCH){

	row row_out;
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

	ALLELE allele_key(allele);
	POLYLL polyll_key(maxll-monoll);
	HWELL hwell_key(maxll-hwell);
	NCUT ncut_key(excluded);
	//NSAMP nsamp_key(site.get); 
	EFCHROM efchrom_key( efc(site) ); 
	MAXLL maxll_key(maxll);

	row_out.add_data((data *)&allele_key);
	row_out.add_data(&polyll_key);
	row_out.add_data(&hwell_key);
	row_out.add_data(&ncut_key);
	row_out.add_data(&efchrom_key);
	row_out.add_data(&maxll_key);
	return row_out;
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

	bool verbose=false;
	bool quite=false;
	bool noheader=false;

	float_t EMLMIN=0.001;
	count_t MIN=0;
	float_t A=0.00;
	float_t MINGOF=2.00;
	count_t MAXPITCH=96;

	count_t skip=0;
	count_t stop=-1;

	std::vector <int> ind;

	/* sets up the help messages and options, see the 'interface.h' for more detials. */

	env_t env;
	env.setname("mapgd ei");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman and Takahiro Maruki");
	env.setdescription("Uses a maximum likelihood approach to estimate population genomic statistics from an individually 'labeled' population.");

	env.optional_arg('i',"input", 	&infile,	&arg_setstr, 	"an error occured while setting the name of the input file.", "the input file for the program (default stdout).");
	env.optional_arg('o',"output", 	&outfile,	&arg_setstr, 	"an error occured while setting the name of the output file.", "the output file for the program (default stdin).");
	env.optional_arg('I',"individuals", &ind, 	&arg_setvectorint, "please provide a list of integers", "the individuals to be used in estimates.\n\t\t\t\ta comma seperated list containing no spaces, and the format X-Y can be used to specify a range (defualt ALL).");
	env.optional_arg('m',"minerror", &EMLMIN, 	&arg_setfloat_t, "please provide a float.", "prior estimate of the error rate (defualt 0.001).");
	env.optional_arg('M',"mincoverage", &MIN, 	&arg_setint, 	"please provide an int.", "minimum coverage of sites to be estimated (defualt 4).");
	env.optional_arg('a',"alpha", 	&A, 		&arg_setfloat_t, "please provide a float.", "cut-off value for printing polymorphic sites (default 0.0).");
	env.optional_arg('g',"goodfit", &MINGOF,	&arg_setfloat_t, "please provide a float.", "cut-off value for the goodness of fit statistic (defaults 2.0).");
	env.optional_arg('N',"number", 	&MAXPITCH,	&arg_setint, 	"please provide an int.", "maximum number of clones to be trimmed (default 96).");
	env.optional_arg('S',"skip", 	&skip,		&arg_setint, 	"please provide an int.", "number of sites to skip before analysis begins (default 0).");
	env.optional_arg('T',"stop", 	&stop,		&arg_setint, 	"please provide an int.", "maximum number of sites to be analyzed (default All sites)");
	env.flag(	'H',"noheader", &noheader,	&flag_set, 	"takes no argument", "disables printing a headerline.");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occured while displaying the version message.", "prints the program version");
	env.flag(	'V', "verbose", &verbose,	&flag_set, 	"an error occured while enabeling verbose excecution.", "prints more information while the command is running.");
	env.flag(	'q', "quite", 	&quite,		&flag_set, 	"an error occured while enabeling quite execution.", "prints less information while the command is running.");

	if ( parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.

	map_file in, out;		//A mapfile is the class that lets us read and write data to a tab delimited file.

					//mapfiles can contain any kind of data, so we must check to see that it has the data 
					//we want, end then we need to format the data for output.

	if (infile.size()!=0) {					//Iff a filename has been set for infile
		in.open(infile.c_str(), std::fstream::in);	
		if (!in.is_open() ) {				//try to open a profile of that name.
			printUsage(env);			//Print help message on failure.
		} 
	} else {
		in.open(std::fstream::in);			//Iff no filename has been set for infile, open profile from stdin.
	};

	if (outfile.size()!=0) {				//Same stuff as infile, now with outfile.	
		out.open(outfile.c_str(), std::fstream::out);
		if (!out.is_open() ){
			printUsage(env);
		}
	} else {
		out.open(std::fstream::out);
	};


	models model;						//model is needed to do our likelihood calculations.

	key locus_key(in.get_key("LOCUS") );			// grabs a key with the name "LOCUS" from the input file.
	key rowid_key(in.get_key("ROWID") );			// grabs a key with the name "ROWID" from the input file.

	std::vector<row> in_row(BUFFER_SIZE, row(in.get_keys() ) );

	std::vector<row> out_row(BUFFER_SIZE, row() );
	

	while (true){			
		uint32_t c=0, readed=0;
		bool estimate_me=1;
		#ifdef PRAGMA
		#pragma omp parallel private(c, model, estimate_me) 
		#endif
		{
			locus this_locus;
			key this_key=locus_key.clone(&this_locus);
			#ifdef PRAGMA
			#pragma omp for
			#endif
			for (uint32_t x=0; x<BUFFER_SIZE; ++x){
				#ifdef PRAGMA
				#pragma omp critical
				#endif
				{
					c=readed;				//Turn on the ability to read data from all clones in 
					if(read_row(in, in_row[c])!=EOF){
						readed++;	//reads the next line of the pro file. pro.read() retuerns 0
						estimate_me=1;
					}
					else estimate_me=0;
				}
				if(estimate_me) {
					get(in_row[c], this_key); 
					out_row[c]=estimate (this_locus, model, MIN, EMLMIN, MINGOF, MAXPITCH);
					
					/*
					#ifdef PRAGMA
					#pragma omp critical
					#endif
					for (size_t i=0; i<sum_gofs.size(); i++){
						sum_gofs[i]+=gofs[i];
						if (gofs[i]!=0) gofs_read[i]++;
					}*/
				}
			}

		}
		real_t bestll, monoll;
		key bestll_key=out_row[0].clone_key("BESTLL", &bestll);
		key monoll_key=out_row[0].clone_key("MONOLL", &monoll);
		if (monoll_key.offset()==key::nokey || bestll_key.offset()==key::nokey) {
			std::cerr << __FILE__ << ":" << __LINE__ << " error: no key returned for row.\n";
			exit(0);
		}
                for (uint32_t x=0; x<readed; ++x){
                	// Now write everything to out. 
			get(out_row[x], bestll_key);
			get(out_row[x], monoll_key);
			if (2*(bestll-monoll)>=A) write_row(out, out_row[x]);
		}
		if (readed!=BUFFER_SIZE){break;}
	}
//	for (size_t x=0; x<ind.size(); ++x)  *out << "@" << pro.get_sample_name(ind[x]) << ":" << sum_gofs[x]/sqrt(float_t(gofs_read[x])) << std::endl;
	in.close();
	out.close();
	return 0;					//Since everything worked, return 0!.
}
