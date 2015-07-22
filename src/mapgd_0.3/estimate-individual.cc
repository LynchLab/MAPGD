/* 

Program estimateIndcpp:

	1) input from a list of site-specific quartets from a labeled population sample;

	2) identify the major and minor alleles, obtain the maximum-likelihood 'like' estimates of allele frequencies, and significance levels.

The designated major nucleotide is simply the one with the highest rank, and the minor nucleotide the one with the second highest rank.
	If the top three ranks are all equal, the site is treated as unresolvable, with both major and minor nucelotides designated in the output by a *.
	If the second and third ranks are equal but lower than the major-nucleotide count, the site is treated as monomorphic, with the minor nucleotide designated by a *.
	
Input File: Three tab delimited dentifiers columns (chromosome, position and ref), followed by an arbitrary number of tab delimited 'quartets', that is counts of the number of times a nucleotide has been observed at a possition.
	Columns are tab delimited, quartets are '/' delimited.
	Default input file is  "datain.txt".

Output File: two columns of site identifiers; reference allele; major allele; minor allele; major-allele frequency; minor-allele frequency; error rate; a ton of other stuff... We really need to clean up the output. 
	Columns are tab delimited.
	Default name is "dataout.txt".
*/

#include "estimate-individual.h"

#define BUFFER_SIZE 5000 
#define PRAGMA

/*
float_t compare (allele_stat mle1, allele_stat mle2, Locus &site1, Locus &site2,  models &model){
	Locus site3=site1+site2;
	alele_stat mle3;
	maximize_grid(site3, mle3, model, gofs, MIN, MAXGOF, MAXPITCH+texc);
	return mle1.ll+mle2.ll-mle3.ll
}*/

/*@breif: Estimates a number of summary statistics from short read sequences.*/ 

//estimate
allele_stat estimate (Locus &site, models &model, std::vector<float_t> &gofs, count_t MIN, count_t EMLMIN, float_t MINGOF, count_t MAXPITCH){

	allele_stat mle, temp;					//allele_stat is a basic structure that containes all the summary statistics for
								//an allele. It gets passed around a lot, and a may turn it into a class that has
								//some basic read and write methods.

	mle.gof=0; mle.efc=0; mle.MM=0; mle.Mm=0; mle.mm=0; 	 //Initialize a bunch of summary statics as 0. 
								 //This should be moved over to the constructor of allele_stat 
								 //(when that constructor is writen). I'm a little concerned that
								 //allele_stat has gotten too bloated, but . . . 
	mle.N=0;

	count_t texc=site.maskedcount();
	
	if (init_params(site, mle, MIN, EMLMIN,0) ){		//If >90% of reads agree, then assume a homozygote,
								//otherwise, assume heterozygote.
	if (mle.error!=0){
		maximize_grid(site, mle, model, gofs, MIN, MINGOF*-1, MAXPITCH+texc);	//trim bad clones and re-fit the model.
//		maximize_newton(site, mle, model, gofs, MIN, MAXGOF, MAXPITCH+texc);	//trim bad clones and re-fit the model.
	}
	else
		maximize_analytical(site, mle, model, gofs, MIN, MINGOF*-1, MAXPITCH+texc);	//trim bad clones and re-fit the model.
	}

	// CALCULATE THE LIKELIHOODS 

	mle.ll=model.loglikelihood(site, mle, MIN);		//Sets the site.ll to the log likelihood of the best fit (ll). 

	if (mle.freq<0.5){					//Check to see if the major and minor alleles are reversed.
		std::swap(mle.major, mle.minor);
		std::swap(mle.MM, mle.mm);
		mle.freq=1.-mle.freq;
		site.swap(0, 1);
		
	}
	else if (mle.freq==0.5){
		if (rand() % 2){				//If the major and minor allele frequencies are identical, 
			std::swap(mle.major, mle.minor);	//flip a coin to determine the major and minor allele.
			std::swap(mle.MM, mle.mm);
			mle.freq=1.-mle.freq;
			site.swap(0, 1);
		};
	};

	temp=mle; temp.MM=1.0; temp.Mm=0.; temp.mm=0.;		//Copies site to mono, then sets mono to a monomophic site 
								//(i.e. sets the genotypic frequencies Mm and mm to 0.
	temp.error=mle.null_error;				//Sets the error rate of mono to the null error rate.
	temp.freq=1.;
	temp.f=0.;
	mle.monoll=model.loglikelihood(site, temp, MIN);			//Sets the site.ll to the log likelihood of the best fit (ll). 

	if (mle.monoll>mle.ll){
		mle.ll=mle.monoll;
		mle.error=mle.null_error;
	};
	temp=mle; 
	temp.MM=pow(mle.freq, 2);				//Similar set up to mono, but now assuming 
	temp.Mm=2.*mle.freq*(1.-mle.freq); 			//Hardy-Weinberg equilibrium.
	temp.mm=pow(1.-mle.freq, 2);				//?
	mle.hwell=model.loglikelihood(site, temp, MIN);		//?
	return mle;
}

int estimateInd(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string infile="";
	std::string outfile="";
	std::string outfilepro;

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
	env.optional_arg('p',"out-pro", &outfilepro,	&arg_setstr, 	"an error occured while setting the name of the output file.", "name of a 'cleaned' pro file (default none).");
	env.optional_arg('I',"individuals", &ind, 	&arg_setvectorint, "please provide a list of integers", "the individuals to be used in estimates.\n\t\t\t\ta comma seperated list containing no spaces, and the format X-Y can be used to specify a range (defualt ALL).");
	env.optional_arg('m',"minerror", &EMLMIN, 	&arg_setfloat_t, "please provide a float.", "prior estimate of the error rate (defualt 0.001).");

//	env.optional_arg('c',"columns", &EMLMIN, 	&arg_setfloat_t, "please provide a float.", "number of columsn in profile (if applicable).");

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

	profile pro, pro_out;		//profile is a fairly complete class that lets us read and write from pro files, 
					//which are files containing set of read 'quartets' that specify the number of 
					//A,C,G and T read at some specific location in a genome. See proFile.h for more info.

	//gcfile out;
	std::ostream *out=&std::cout;
	std::ofstream outFile;

	genotype dip(2);
	genotype trip(3);
	genotype tet(4);

	if (infile.size()!=0) {					//Iff a filename has been set for infile
		if (pro.open(infile.c_str(), "r")==NULL) {	//try to open a profile of that name.
			printUsage(env);			//Print help message on failure.
		} 
	}
	else {
		pro.open("r");					//Iff no filename has been set for infile, open profile from stdin.
	};

	if (outfile.size()!=0) {
		outFile.open(outfile, std::ofstream::out);
		if (!outFile) printUsage(env);
		out=&outFile;
	};

	//else out.open('w', CSV);				//Iff no filename has been set for outfile, pgdfile prints to stdout.

	count_t outc=6;
	char cdel='\t';
	char qdel='/';
	bool binary=false;
	if (outfilepro.size()!=0) {				//Same sort of stuff for the outf
		if (binary) {
			pro_out.open(outfilepro.c_str(), "wb");
			if (!pro_out.is_open() ) {printUsage(env); exit(0);}
		} else {
			pro_out.open(outfilepro.c_str(), "w");
			if (!pro_out.is_open() ) {printUsage(env); exit(0);}
		}
		pro_out.setsamples(pro.size() );
		pro_out.setcolumns(outc);
		pro_out.set_delim_column(cdel);
		pro_out.set_delim_quartet(qdel);
		if (not (noheader) ) pro_out.writeheader();
	};

	/* this is the basic header of our outfile, should probably be moved over to a method in allele_stat.*/
	if (not (noheader) ) *out << "id1\tid2\tref\tmajor\tminor\tcov\tM\tm\terror\tnull_e\tf\tMM\tMm\tmm\th\tpolyll\tHWEll\tgof\teff_chrom\tN\tN_excluded\tmodel_ll" << std::endl;

	pro.maskall();							//Turn off the ability to read data from all clones by default. 

	if ( ind.size()==0 ) { 						//Iff the vector ind (which should list the clones to 
		ind.clear();						//be read from the .pro file) is empty, then 
		for (count_t x=0; x<pro.size(); ++x) ind.push_back(x);  //put every clone in the vector ind.
	};

	std::vector <float_t> sum_gofs(ind.size() );
	std::vector <float_t> gofs_read(ind.size() );
	models model;
	allele_stat buffer_mle[BUFFER_SIZE]; 
	Locus buffer_site[BUFFER_SIZE];
	while (true){			//reads the next line of the pro file. pro.read() retuerns 0
		uint32_t c=0, readed=0;
		bool estimate_me=1;
		#ifdef PRAGMA
		#pragma omp parallel private(c, model, estimate_me) 
		#endif
		{
			#ifdef PRAGMA
			#pragma omp for
			#endif
			for (uint32_t x=0; x<BUFFER_SIZE; ++x){
				#ifdef PRAGMA
				#pragma omp critical
				#endif
				{
					c=readed;				//Turn on the ability to read data from all clones in 
					if(pro.read(buffer_site[c])!=EOF){
						readed++;	//reads the next line of the pro file. pro.read() retuerns 0
						estimate_me=1;
					}
					else estimate_me=0;
				}
				if(estimate_me) {
					std::vector <float_t> gofs(ind.size() );
					buffer_mle[c]=estimate (buffer_site[c], model, gofs, MIN, EMLMIN, MINGOF, MAXPITCH);
					#ifdef PRAGMA
					#pragma omp critical
					#endif
					for (size_t i=0; i<sum_gofs.size(); i++){
						sum_gofs[i]+=gofs[i];
						if (gofs[i]!=0) gofs_read[i]++;
					}
				}
			}

		}
                for (uint32_t x=0; x<readed; ++x){
                	// Now print everything to the *out stream, which could be a file or the stdout. 
			//TODO move this over into a formated file.
			//?
			if (2*(buffer_mle[x].ll-buffer_mle[x].monoll)>=A){
				*out << std::fixed << std::setprecision(6) << pro.getids(buffer_site[x]) << '\t' << buffer_site[x].getname(0) << '\t' << buffer_site[x].getname_gt(1) << '\t';
				*out << std::fixed << std::setprecision(6) << buffer_mle[x] << std::endl;
			}
			if (buffer_mle[x].gof<-MINGOF) buffer_site[x].maskall(); 
			if (pro_out.is_open() ){
				buffer_site[x].id0=pro_out.encodeid0(pro.decodeid0(buffer_site[x].id0) );
				pro_out.write(buffer_site[x]);
			}
		}
		if (readed!=BUFFER_SIZE){break;}
	}
	for (size_t x=0; x<ind.size(); ++x)  *out << "@" << pro.getsample_name(x) << ":" << sum_gofs[x]/sqrt(float_t(gofs_read[x])) << std::endl;
	pro.close();
	if (outFile.is_open()) outFile.close();		//Closes outFile iff outFile is open.
	if (pro_out.is_open()) pro_out.close();		//Closes pro_out iff pro_out is open.
	return 0;					//Since everything worked, return 0!.
}
