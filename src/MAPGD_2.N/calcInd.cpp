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

#include "calcInd.h"

/*@breif: Estimates a number of summary statistics from short read sequences.*/ 

int calcInd(int argc, char *argv[])
{
	std::string infile="datain.txt";

	/* All the variables that can be set from the command line */

	float_t EMLMIN=0.001;
	count_t MIN=4;
	float_t a=0.00, MM=1.0, Mm=0.0, mm=0;
	float_t maxgof=2.00;
	count_t maxpitch=96;

	std::vector <int> ind;
	std::vector <std::string> loc;

	/* sets up the help messages and options, see the 'interface.h' for more detials. */

	env_t env;
	env.setname("mapgd calc");
	env.setver("0.0");
	env.setauthor("Matthew Ackerman");
	env.setdescription("Calculates the log likelihood of a line in a pro file. This command is currently unstable.");

	//env.required_arg('g',"goto", &loc, 	&arg_setvectorstr, "please provide a colon seperated locatation", "The location in the pro file at which to caclulate the ll stat.");
	env.required_arg('P',"individuals", &MM, 	&arg_setfloat_t, "please provide a list of integers", "Choose individuals to use in estimate (defualts to ALL, currently broken).");
	env.required_arg('H',"individuals", &Mm, 	&arg_setfloat_t, "please provide a list of integers", "Choose individuals to use in estimate (defualts to ALL, currently broken).");


	env.optional_arg('I',"individuals", &ind, 	&arg_setvectorint, "please provide a list of integers", "Choose individuals to use in estimate (defualts to ALL, currently broken).");
	env.optional_arg('e',"error", &EMLMIN, 	&arg_setfloat_t, "please provide a float.", "prior estimate of the error rate (defualt 0.001).");
	env.optional_arg('M',"mincoverage", &MIN, 	&arg_setint, 	"please provide an int.", "minimum coverage of sites to be estimated (defualt 4).");
	env.optional_arg('a',"alpha", 	&a, 		&arg_setfloat_t, "please provide a float.", "cut-off value for printing polymorphic sites (default 0.0).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occured while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occured while displaying the version message.", "prints the program version");

	if ( parsargs(argc, argv, env) ) printUsage(env); //Gets all the command line options, and prints usage on failure.

	std::ostream *out;		//declare the stream will will write to.
	out=&std::cout;			//set the outstream to the std:cout by default.

	profile pro;			//profile is a fairly complete class that lets us read and write from pro files, 
					//which are files containing set of read 'quartets' that specify the number of 
					//A,C,G and T read at some specific location in a genome. See proFile.h for more info.

	allele_stat mle;		//allele_stat is a basic structure that containes all the summary statistics for
					//an allele. It gets passed around a lot, and a may turn it into a class that has
					//some basic read and write methods.

	mle.error=EMLMIN; mle.efc=0; mle.MM=MM; mle.Mm=Mm; mle.mm=1.-MM-mm; //Initialize a bunch of summary statics as 0. 
								 //This should be moved over to the constructor of allele_stat 
								 //(when that constructor is writen). I'm a little concerned that
								 //allele_stat has gotten too bloated, but . . . 

	
	if (infile.size()!=0) {					//Iff a filename has been set for infile
		if (pro.open(infile.c_str(), 'r')==NULL) {	//try to open a profile of that name.
			printUsage(env);			//Print help message on failure.
		} 
	}
	else pro.open('r');					//Iff no filename has been set for infile, open profile from stdin.

	/* this is the basic header of our outfile, should probably be moved over to a method in allele_stat.*/
	*out << "#id1\tid2\tref\tmajor\tminor\tpopN\tM\tm\terror\tnull_error\th\tpolyll\tHWEll\tgof\tefN\tN_excluded" << std::endl;

	pro.maskall();	//Turn off the ability to read data from all clones by default. 

	if ( ind.size()==0 ) { 						//Iff the vector ind (which should list the clones to 
		ind.clear();						//be read from the .pro file) is empty, then 
		for (count_t x=0; x<pro.size(); ++x) ind.push_back(x);  //put every clone in the vector ind.
	};

	for (count_t x=0; x<ind.size(); ++x) pro.unmask(ind[x]);	//Turn on the ability to read data from all clones in 
									//the vector ind.

	//pro.seek(loc);
	pro.sort();

	mle.ll=ll(pro, mle, MIN);				//Sets the site.ll to the log likelihood of the best fit (ll). 
		
        /* Now print everything to the *out stream, which could be a file or the stdout. */
	//if (llPOLY>a){
	*out << std::fixed << std::setprecision(6) << pro.getids() << '\t' << pro.getname(0) << '\t' << pro.getname_gt(1) << '\t' << '\t' << mle.ll << std::endl << '\t';
	pro.close();					//Close the pro file/stream.
	exit(0);					//Since everything worked, return 0!.
}
