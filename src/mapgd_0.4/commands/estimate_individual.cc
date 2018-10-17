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

#include "estimate_individual.h"

#define BUFFER_SIZE 500

/*
float_t compare (Allele mle1, Allele mle2, Locus &site1, Locus &site2,  models &model){
	Locus site3=site1+site2;
	alele_stat mle3;
	maximize_grid(site3, mle3, model, gofs, MIN, MAXGOF, MAXPITCH+texc);
	return mle1.ll+mle2.ll-mle3.ll
}*/

/// Estimates a number of summary statistics from short read sequences.
/**
 */
Allele estimate (Locus &site, models &model, std::vector<float_t> &gofs, const count_t &MIN, const float_t &EMLMIN, const float_t &MINGOF, const count_t &MAXPITCH, bool newton, bool bias){

#ifdef MPI
	std::cerr << "HOW THE HELL IS THIS HAPPENING!\n";
	exit(0);
	MPI_Init(&argc, &argv);

	int numtasks, taskid;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
#endif


	Allele mle, temp;					//Allele is a basic structure that containes all the summary statistics for
								//an allele. It gets passed around a lot, and a may turn it into a class that has
								//some basic read and write methods.

	mle.gof=0; mle.efc=0; mle.MM=0; mle.Mm=0; mle.mm=0; 	 //Initialize a bunch of summary statics as 0. 
								 //This should be moved over to the constructor of Allele 
								 //(when that constructor is writen). I'm a little concerned that
								 //Allele has gotten too bloated, but . . . 
	mle.N=0;
	mle.set_abs_pos(site.get_abs_pos() );
	mle.ref=site.ref.base;
	mle.major=4;
	mle.minor=4;

	site.mask_low_cov(MIN-1);
	count_t texc=site.maskedcount(), rexc;
	rexc=texc;

	if (init_params(site, mle, EMLMIN) ){		//If >90% of reads agree, then assume a homozygote,
							//otherwise, assume heterozygote.
	if (mle.null_error!=0){
		rexc=maximize_grid(site, mle, model, gofs, -MINGOF, MAXPITCH+texc);	//trim bad clones and re-fit the model.
		if (newton) 
			rexc=maximize_restricted_newton(site, mle, model, gofs, -MINGOF, MAXPITCH+texc);	//the NR maximization
	}
	else
		rexc=maximize_analytical(site, mle, model, gofs, -MINGOF, MAXPITCH+texc);//trim bad clones and re-fit the model.
	} else {
		return mle;
	}

	mle.excluded=rexc-texc;

	// CALCULATE THE LIKELIHOODS 

	mle.ll=model.loglikelihood(site, mle);			//Sets the site.ll to the log likelihood of the best fit (ll). 


	if (mle.freq<0.5){					//Check to see if the major and minor alleles are reversed.
		std::swap(mle.major, mle.minor);
		std::swap(mle.MM, mle.mm);
		std::swap(mle.null_error, mle.null_error2);
		mle.freq=1.-mle.freq;
		site.swap(0, 1);
		
	} 
	/*else if (mle.freq==0.5){
		//struct drand48_data drand_buf;
		//drand48_r (&drand_buf, &r);
		if (r % 2){					//If the major and minor allele frequencies are identical, 
			std::swap(mle.major, mle.minor);	//flip a coin to determine the major and minor allele.
			std::swap(mle.MM, mle.mm);
			std::swap(mle.null_error, mle.null_error2);
			mle.freq=1.-mle.freq;
			site.swap(0, 1);
		};
	};*/

	temp=mle; temp.MM=1.0; temp.Mm=0.; temp.mm=0.;		//Copies site to mono, then sets mono to a monomophic site 
								//(i.e. sets the genotypic frequencies Mm and mm to 0.
	temp.error=mle.null_error;				//Sets the error rate of mono to the null error rate.
	temp.freq=1.;
	temp.f=0.;

	float_t mono1, mono2;

	mono1=model.loglikelihood(site, temp);			//Sets the site.ll to the log likelihood of the best fit (ll). 
	
	temp=mle; temp.MM=0.0; temp.Mm=0.; temp.mm=1.;		
	temp.error=mle.null_error2;				//Sets the error rate of mono to the null error rate.
	temp.freq=0.;
	temp.f=0.;

	mono2=model.loglikelihood(site, temp);			//Sets the site.ll to the log likelihood of the best fit (ll). 
	mono1 > mono2 ? mle.monoll=mono1 : mle.monoll=mono2; 

	if (mle.monoll>mle.ll){
		mle.ll=mle.monoll;
		mle.error=mle.null_error;
	};

	temp=mle; 
	temp.MM=pow(mle.freq, 2);				//Similar set up to mono, but now assuming 
	temp.Mm=2.*mle.freq*(1.-mle.freq); 			//Hardy-Weinberg equilibrium.
	temp.mm=pow(1.-mle.freq, 2);				//?
	mle.hwell=model.loglikelihood(site, temp);		//?
	
	if (bias){
		get_bias(site, mle);
		mle.print_bias=true;
	}

	/*if (site.getcount(0)==site.getcount(1) ){
		mle.major=4;
		mle.minor=4;
	} else */
	if (site.getcount(1)==site.getcount(2) ) mle.minor=4;
	return mle;
}

#ifdef MPI
void write(std::ostream& out, uint32_t readed, profile& pro, profile& pro_out, Allele* buffer_mle, Locus* buffer_site, float MINGOF, float_t A)
{
	for (uint32_t x = 0; x < readed; ++x) {
		// Now print everything to the *out stream, which could be a file or the stdout. 
		//TODO move this over into a formated file.
		//?
		if (2 * (buffer_mle[x].ll - buffer_mle[x].monoll) >= A) {
			out << std::fixed << std::setprecision(6) << pro.getids(buffer_site[x]) << '\t' << buffer_site[x].getname(0) << '\t' << buffer_site[x].getname_gt(1) << '\t';
			out << std::fixed << std::setprecision(6) << buffer_mle[x] << std::endl;
		}
		if (buffer_mle[x].gof < -MINGOF) buffer_site[x].maskall();
		if (pro_out.is_open()) {
			buffer_site[x].id0 = pro_out.encodeid0(pro.decodeid0(buffer_site[x].id0));
			pro_out.write(buffer_site[x]);
		}
	}
}

bool mpi_select(int lineid, int taskid, int num_tasks)
{
	lineid = lineid % (BUFFER_SIZE*num_tasks);
	return lineid >= taskid*BUFFER_SIZE && lineid < (taskid + 1)*BUFFER_SIZE;
}

std::string mpi_recieve_string(int rank)
{
	MPI_Status status;
	MPI_Probe(rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	int count;
	MPI_Get_count(&status, MPI_CHAR, &count);
	char *buf = new char[count];
	MPI_Recv(buf, count, MPI_CHAR, rank, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	std::string result(buf, count);
	delete[] buf;
	return result;
}

void 
do_estimate(Allele* buffer_mle, Locus& buffer_site, models& model, 
	std::vector<int>& ind, std::vector <float_t>& sum_gofs, 
	std::vector <float_t>& gofs_read, count_t MIN, float_t EMLMIN, 
	float_t MINGOF, count_t MAXPITCH)
{
	std::vector <float_t> gofs(ind.size());
	buffer_site.unmaskall();
	*buffer_mle = estimate(buffer_site, model, gofs, MIN, EMLMIN, MINGOF, MAXPITCH);
	if (2 * (buffer_mle->ll - buffer_mle->monoll) >= 22) {
		for (size_t i = 0; i < sum_gofs.size(); i++) {
			sum_gofs[i] += gofs[i];
			if (gofs[i] != 0) gofs_read[i]++;
		}
	}

}
#endif

int estimateInd(int argc, char *argv[])
{

	/* All the variables that can be set from the command line */

	std::string infile="";
	std::string outfile="";
	std::string outfilepro;
	std::string indexname="";

	bool verbose=false;
	bool quite=false;
	bool noheader=false;
	bool newton=false;
	bool bias=false;

	int rnseed=3;

	double EMLMIN=0.0001;
	int MIN=4;
	double MINGOF=2.00;
	int MAXPITCH=96;

	std::vector <size_t> ind;

	/* sets up the help messages and options, see the 'interface.h' for more detials. */
	//std::cerr "If you are using .... please cite ...."

	Environment env;
	env.set_name("mapgd allele");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman and Takahiro Maruki");
	env.set_description("Uses a maximum likelihood approach to estimate population genomic statistics from an individually 'labeled' population.");
	env.optional_arg('r',"seed", 	rnseed,	"please provide a number.", "random number seed (3).");

	env.optional_arg('o',"output", 	outfile,	"an error occurred while setting the name of the output file.", "the output file for the program (default stdin).");
	env.optional_arg('p',"outpro",  outfilepro,	"an error occurred while setting the name of the output file.", "name of a 'cleaned' pro file (default none).");
	env.optional_arg('I',"individuals", ind, 	"please provide a list of integers", "the individuals to be used in estimates. A comma separated list containing no spaces, the python slice notation (i.e. 1:4 for 1,2,3,4 or [1:5:2] for 1,3,4 ... can be used to specify this list (default ALL).");
	env.optional_arg('e',"min-error", EMLMIN, 	"please provide a float.", "prior estimate of the error rate (default 0.001).");

//	env.optional_arg('c',"columns", &EMLMIN, 	&arg_setfloat_t, "please provide a float.", "number of columsn in profile (if applicable).");

	env.optional_arg('H',"header", indexname, 	"please provide an str.", "the name of a .idx file storing scaffold information");
	env.optional_arg('c',"min-coverage", MIN, 	"please provide an int.", "minimum coverage for an individual at a site for an individual to be used (default 4).");
	env.optional_arg('g',"good-fit", MINGOF,	"please provide a float.", "cut-off value for the goodness of fit statistic (defaults 2.0).");
	env.optional_arg('B',"max-bad",  MAXPITCH,	"please provide an int.", "cut-off value for number of bad individuals needed before a site is removed entirely (default 96).");

	env.positional_arg('i',"input",	infile,	"No input file specified", "the input file for the program (default stdout).");

	env.flag(	'b',"bias",  &bias,	&flag_set, 	"takes no argument", "print major allele bias.");
	env.flag(	'N',"noheader", &noheader,	&flag_set, 	"takes no argument", "disables printing a header-line.");
	env.flag(	'n',"newton", 	&newton,	&flag_set, 	"takes no argument", "use Newton-Raphson likelihood maximization (not working).");
	env.flag(	'h', "help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message.", "prints this message");
	env.flag(	'v', "version", &env, 		&flag_version, 	"an error occurred while displaying the version message.", "prints the program version");
	env.flag(	'V', "verbose", &verbose,	&flag_set, 	"an error occurred while enabling verbose execution.", "prints more information while the command is running.");
	env.flag(	'q', "quiet", 	&quite,		&flag_set, 	"an error occurred while enabling quite execution.", "prints less information while the command is running.");


	if ( parsargs(argc, argv, env) ) print_usage(env); //Gets all the command line options, and prints usage on failure.

	Indexed_file <Allele> map_out;
	Indexed_file <Locus> pro_in, pro_out;
	
	Allele allele_out;
	Locus locus_out, locus_in;

	if (infile.size()!=0) {					//Iff a filename has been set for infile
		pro_in.open(infile.c_str(), std::fstream::in);	
		if (!pro_in.is_open() ) {				//try to open a profile of that name.
			print_usage(env);			//Print help message on failure.
		} 
		locus_in=pro_in.read_header();
	}
	else {
		pro_in.open(std::fstream::in);			//Iff no filename has been set for infile, open profile from stdin.
		locus_in=pro_in.read_header();
	};
	if (!pro_in.table_is_open() ) 
	{
		fprintf(stderr, "%s:%d. Error: cannot open file.\n",__FILE__, __LINE__);
		exit(0);
	}

	bool binary=false;

	if (outfilepro.size()!=0) {				//Same sort of stuff for the outf
		if (binary) {
			pro_out.open(outfilepro.c_str(), std::fstream::out | std::fstream::binary);
			if (!pro_out.is_open() ) {print_usage(env); exit(0);}
		} else {
			pro_out.open(outfilepro.c_str(), std::fstream::out);
			if (!pro_out.is_open() ) {print_usage(env); exit(0);}
		}
		locus_out.set_sample_names(locus_in.get_sample_names() );
		pro_out.set_index(pro_in.get_index() );
		pro_out.write_header(locus_out);
	};


	/* this is the basic header of our outfile, should probably be moved over to a method in Allele.*/
	locus_in.maskall();						//Turn off the ability to read data from all clones by default. 

	if ( ind.size()==0 ) { 						//Iff the vector ind (which should list the clones to 
		ind.clear();						//be read from the .pro file) is empty, then 
		for (count_t x=0; x<locus_in.get_sample_names().size(); ++x) ind.push_back(x);  //put every clone in the vector ind.
	};

	std::vector <float_t> sum_gofs(ind.size() );
	std::vector <float_t> gofs_read(ind.size() );

	models model;

	if (outfile.size()!=0) {
		map_out.open(outfile.c_str(), std::ios::out);
		if (!map_out.is_open() ) print_usage(env);
	} else 	map_out.open(std::fstream::out);

	Allele buffer_mle[BUFFER_SIZE]; 
	Locus buffer_locus[BUFFER_SIZE];
	std::fill_n(buffer_locus, BUFFER_SIZE, locus_in);

	map_out.set_index(pro_in.get_index() );
	map_out.write_header(buffer_mle[0]);

	uint32_t all_read=0;

	while (true){			//reads the next line of the pro file. pro.read() retuerns 0
		uint32_t c=0, readed=0;
		bool estimate_me=true;
		#ifndef NOOMP
		#pragma omp parallel private(c, model, estimate_me) 
		#endif
		{
			#ifndef NOOMP
			#pragma omp for
			#endif
			for (uint32_t x=0; x<BUFFER_SIZE; ++x){
				#ifndef NOOMP
				#pragma omp critical
				#endif
				{
					c=readed;		//Turn on the ability to read data from all clones in 
					if(pro_in.read(buffer_locus[c]).table_is_open() ){
						readed++;	//reads the next line of the pro file. pro.read() retuerns 0
						estimate_me=1;
					}
					else estimate_me=0;
				}
				if(estimate_me) {
					std::vector <float_t> gofs(ind.size() );
		//			buffer_locus[c].unmaskall();
					buffer_locus[c].unmask(ind);
					buffer_mle[c]=estimate(buffer_locus[c], model, gofs, MIN, EMLMIN, MINGOF, MAXPITCH, newton, bias);

					#ifndef NOOMP
					#pragma omp critical
					#endif
					if (2*(buffer_mle[c].ll-buffer_mle[c].monoll)>=22 && buffer_mle[c].gof>=-2 ){
						for (size_t i=0; i<sum_gofs.size(); i++){
							sum_gofs[i]+=gofs[i];
							if (gofs[i]!=0) gofs_read[i]++;
						}
					}
				}
			}

		}
                for (uint32_t x=0; x<readed; ++x){
			map_out.write(buffer_mle[x]);
			if (buffer_mle[x].gof<-MINGOF) buffer_locus[x].maskall(); 
			if (pro_out.is_open() ){
				pro_out.write(buffer_locus[x]);
			}
		}
		if (readed!=BUFFER_SIZE){break;}
		all_read+=readed;
	}
	map_out.close_table();
	if ( not(noheader) ) {
		Flat_file <Sample_gof> gof_file;
		gof_file.open_from(map_out);
		gof_file.write_header(Sample_gof() );
		for (size_t x=0; x<ind.size(); ++x) gof_file.write(Sample_gof(locus_in.get_sample_names()[ind[x]], sum_gofs[x]/(float_t(gofs_read[x]) ) ) );
		gof_file.close();
	}
	pro_in.close();
	if (pro_out.is_open()) pro_out.close();		//Closes pro_out iff pro_out is open.
	return 0;					//Since everything worked, return 0!.
}
