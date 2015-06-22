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

#include "estimateInd.h"

/*@breif: Estimates a number of summary statistics from short read sequences.*/ 

//estimate
//print

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
	count_t MIN=4;
	float_t a=0.00;
	float_t maxgof=2.00;
	float_t	tgof=0;
	count_t maxpitch=96;

	count_t skip=0;
	count_t stop=-1;

	std::vector <int> ind;

	/* sets up the help messages and options, see the 'interface.h' for more detials. */
	infile=""; outfile=""; 

	env_t env;
	env.setname("mapgd ei");
	env.setver("2.1");
	env.setauthor("Takahiro Maruki and Matthew Ackerman");
	env.setdescription("Uses a maximum likelihood approach to estimate population genomic statistics from an individually 'labeled' population.");

	env.optional_arg('i',"input", 	&infile,	&arg_setstr, 	"an error occured while setting the name of the input file.", "sets the input file for the program (default 'datain.txt').");
	env.optional_arg('o',"output", 	&outfile,	&arg_setstr, 	"an error occured while setting the name of the output file.", "sets the output file for the program (default 'dataout.txt').");
	env.optional_arg('p',"out-pro", &outfilepro,	&arg_setstr, 	"an error occured while setting the name of the output file.", "sets the output file for the program (default 'dataout.txt').");
	env.optional_arg('I',"individuals", &ind, 	&arg_setvectorint, "please provide a list of integers", "Choose individuals to use in estimate.\n\t\t\t\t Should be a comma seperated list containing no spaces, and the format X-Y can be used to specify a range (defualts to ALL).");
	env.optional_arg('m',"minerror", &EMLMIN, 	&arg_setfloat_t, "please provide a float.", "prior estimate of the error rate (defualt 0.001).");
	env.optional_arg('M',"mincoverage", &MIN, 	&arg_setint, 	"please provide an int.", "minimum coverage of sites to be estimated (defualt 4).");
	env.optional_arg('a',"alpha", 	&a, 		&arg_setfloat_t, "please provide a float.", "cut-off value for printing polymorphic sites (default 0.0).");
	env.optional_arg('g',"goodfit", &maxgof,	&arg_setfloat_t, "please provide a float.", "cut-off value for the goodness of fit statistic (defaults 2.0).");
	env.optional_arg('N',"number", 	&maxpitch,	&arg_setint, 	"please provide an int.", "maximum number of clones to be trimmed (default 96).");
	env.optional_arg('S',"skip", 	&skip,		&arg_setint, 	"please provide an int.", "Number of sites to skip before analysis begins (default 0).");
	env.optional_arg('T',"stop", 	&stop,		&arg_setint, 	"please provide an int.", "Maximum number of sites to be analyzed (default All sites)");
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

	allele_stat mle, hwe, mono;	//allele_stat is a basic structure that containes all the summary statistics for
					//an allele. It gets passed around a lot, and a may turn it into a class that has
					//some basic read and write methods.

	mle.gof=0; mle.efc=0; mle.MM=0; mle.Mm=0; mle.mm=0; 	 //Initialize a bunch of summary statics as 0. 
								 //This should be moved over to the constructor of allele_stat 
								 //(when that constructor is writen). I'm a little concerned that
								 //allele_stat has gotten too bloated, but . . . 
	int line=0;

	if (infile.size()!=0) {					//Iff a filename has been set for infile
		if (pro.open(infile.c_str(), "r")==NULL) {	//try to open a profile of that name.
			printUsage(env);			//Print help message on failure.
		} 
	}
	else pro.open("r");					//Iff no filename has been set for infile, open profile from stdin.

	if (outfile.size()!=0) {
		outFile.open(outfile, std::ofstream::out);
		if (!outFile) printUsage(env);
		out=&outFile;
	};

	//else out.open('w', CSV);				//Iff no filename has been set for outfile, pgdfile prints to stdout.

	//Set up file?
	//out.setheader(); << "id1\tid2\tref\tmajor\tminor\tcov\tM\tm\terror\tnull_e\tf\tMM\tMm\tmm\th\tpolyll\tHWEll\tgof\teff_chrom\tN\tN_excluded\tmodel_ll" << std::endl;

	if (outfilepro.size()!=0) {				//Same sort of stuff for the outfile. 
		pro_out.open(outfilepro.c_str(), "w");
		if (!pro_out.is_open()){
			std::cerr << "Cannot open file " << outfilepro << ". This file may already exist, or you may be trying to write to ro location.";
			printUsage(env);
		};
		pro_out.copyheader(pro);
		pro_out.writeheader();
	};


	/* this is the basic header of our outfile, should probably be moved over to a method in allele_stat.*/
	if (not (noheader) ) *out << "id1\tid2\tref\tmajor\tminor\tcov\tM\tm\terror\tnull_e\tf\tMM\tMm\tmm\th\tpolyll\tHWEll\tgof\teff_chrom\tN\tN_excluded\tmodel_ll" << std::endl;

	pro.maskall();	//Turn off the ability to read data from all clones by default. 

	if ( ind.size()==0 ) { 						//Iff the vector ind (which should list the clones to 
		ind.clear();						//be read from the .pro file) is empty, then 
		for (count_t x=0; x<pro.size(); ++x) ind.push_back(x);  //put every clone in the vector ind.
	};

	for (count_t x=0; x<ind.size(); ++x){
		 pro.unmask(ind[x]);					//Turn on the ability to read data from all clones in 
	}								//the vector ind.

	count_t excluded, texc, read=0;					//A varaiable for keeping track of the number of clones
									//we have excluded from the analysis.
	for (int x=0; x<skip; ++x) pro.read(SKIP);
	

	std::chrono::time_point <std::chrono::system_clock> then, now;
	std::chrono::duration <double> runtime;
	then=std::chrono::system_clock::now();				//Right now it is then, but it will be now later. . .
	float_t secleft;

	count_t nextprint=95;
	while (pro.read()!=EOF ){					//reads the next line of the pro file. pro.read() retuerns 0
									//on success, EOF when end of file reached. (switch to ==0?)

		if ( pro.size()!=0 && outFile.is_open() && !quite && read>nextprint){
			now=std::chrono::system_clock::now();
			runtime=now-then;
			secleft=read/runtime.count();
			//I should probably use std::put_time here, but it is not implemented in some recent versions of GCC, so I will hold off 
			//until later. 
  			std::cout << "Lines per second: " << secleft << std::endl;
			nextprint=read+1000;
		}

		//estimate(line from profile, );
		//print();
		
		mle.N=0;

		texc=pro.maskedcount();

		initparams(pro, mle, MIN, EMLMIN,0);
									//If >90% of reads agree, then assume a homozygote,
									//otherwise, assume heterozygote.

		excluded=maximizegrid(pro, mle, MIN, maxgof, maxpitch+texc);	//trim so clones and re-fit the model.

		tgof=mle.gof; 
									
		excluded=pro.maskedcount();				//count the number of clones excluded. 

                // CALCULATE THE LIKELIHOODS 

		mle.ll=loglikelihood(pro, mle, MIN);			//Sets the site.ll to the log likelihood of the best fit (ll). 
		mono=mle; mono.MM=1.0; mono.Mm=0.; mono.mm=0.;		//Copies site to mono, then sets mono to a monomophic site 
									//(i.e. sets the genotypic frequencies Mm and mm to 0.
		mono.error=mle.null_error;				//Sets the error rate of mono to the null error rate.

		mono.freq=1.;
		mono.f=0.;
                mono.ll=loglikelihood(pro, mono, MIN);				//Calculates the log likelihood of the mono fit.
		if (mono.ll>mle.ll){
			mle=mono;
			maximizegrid(pro, mle, MIN, maxgof, maxpitch+texc);
		};
		mono.ll=(mle.ll-mono.ll)*2.;
	
		if (mle.freq<0.5){					//Check to see if the major and minor alleles are reversed.
			std::swap(mle.major, mle.minor);
			std::swap(mle.MM, mle.mm);
			mle.freq=1.-mle.freq;
			pro.swap(0, 1);
			
		}
		else if (mle.freq==0.5){
			if (rand() % 2){				//If the major and minor allele frequencies are identical, 
				std::swap(mle.major, mle.minor);	//flip a coin to determine the major and minor allele.
				std::swap(mle.MM, mle.mm);
				mle.freq=1.-mle.freq;
				pro.swap(0, 1);
			};
		};

		hwe=mle; 
		hwe.MM=pow(mle.freq, 2);				//Similar set up to mono, but now assuming 
		hwe.Mm=2.*mle.freq*(1.-mle.freq); 			//Hardy-Weinberg equilibrium.
		hwe.mm=pow(1.-mle.freq, 2);				//?
		hwe.ll=(mle.ll-loglikelihood(pro, hwe, MIN))*2;		//?
		
                // Now print everything to the *out stream, which could be a file or the stdout. 
		//TODO move this over into a formated file.
		if (mono.ll>=a){
			if (mle.N>0){
				*out << std::fixed << std::setprecision(6) << pro.getids() << '\t' << pro.getname(0) << '\t' << pro.getname_gt(1) << '\t' << pro.getcoverage() << '\t' << mle.freq <<'\t' << 1.-mle.freq << '\t' << mle.error << '\t';
				*out << std::fixed << std::setprecision(6) << mle.null_error <<'\t' << mle.f << '\t' << mle.MM << '\t' << mle.Mm <<'\t' << mle.mm << '\t' << mle.freq*(1.-mle.freq)*2 << '\t' << mono.ll << '\t' << hwe.ll << '\t' << tgof << '\t' << mle.efc << '\t' << mle.N << '\t' << excluded-texc << '\t' << -2*mle.ll << '\n';
			}
			else{
				*out << std::fixed << std::setprecision(6) << pro.getids() << '\t' << '*' << '\t' << '*' << '\t' << pro.getcoverage() << '\t' << '*' <<'\t' << '*' << '\t' << '*' << '\t';
				*out << std::fixed << std::setprecision(6) << '*' <<'\t' << '*' << '\t' << '*' << '\t' << '*' <<'\t' << '*' << '\t' << '*' << '\t' << '*' << '\t' << '*' << '\t' << '*' << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\t' << '*' << '\n';
			}
		}

		if (pro_out.is_open() ){
			pro_out.copy(pro);
			if (tgof<-maxgof) pro_out.maskall(); 
			pro_out.write();
		};

		for (count_t x=0; x<ind.size(); ++x) pro.unmask(ind[x]);	//Turn on the ability to read data from all clones in 
		if (read==stop) break;
		read++;
	}
	pro.close();					//Close the pro file/stream.
	if (outFile.is_open()) outFile.close();		//Closes outFile iff outFile is open.
	if (pro_out.is_open()) pro_out.close();		//Closes pro_out iff pro_out is open.
	exit(0);					//Since everything worked, return 0!.
}
