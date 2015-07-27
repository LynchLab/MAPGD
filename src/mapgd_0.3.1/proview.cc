/******************* proview.cpp *******************
 * This is inspired by sam2pro, some of the libraries
 * by Bernhard Haubold, haubold@evolbio.mpg.de
 * 			and
 * by Matthew Ackerman, matthew.s.ackerman@gmail.com
 * Description: Convert sam output to profiles.
 * Date: Wed Jul 17 12:00:00 2015
 **************************************************/

//TODO: change merging behavior to merge samples.

#include "proview.h"

int proview(int argc, char *argv[])
{
	/* commands that can be set from the command line */

	/* default values for arguments */
	Args args;

	char cdel='\t';
	char qdel='/';

	bool noheader=false;
	bool binary=false;
	bool bernhard=false;
	bool split=false;

	count_t outc=6;
	count_t inC=0;

	args.pro=false;
	args.min=4;
	args.pvalue=0.001;

	std::vector <std::string> infiles;
	int sample=-1;
	std::string outfile="";

	env_t env;
	env.setname("mapgd proview");
	env.setver(VERSION);
	env.setauthor("Matthew Ackerman and Bernhard Haubold");
	env.setdescription("prints data in the '.pro' file quartet format");
	
	env.optional_arg('m',"minimum",	 &args.min, 	&arg_setint, 	"an error occured", "prints a line iff at least one line has coverage greater than the minimum coverage (defauld 4)");
	env.optional_arg('q',"qdel",	&qdel,		&arg_setchar, 	"an error occured", "sets the output quartet delimiter (default tab)");
	env.optional_arg('d',"cdel",	&cdel,		&arg_setchar, 	"an error occured", "sets the output column delimiter (default tab)");
	env.optional_arg('c',"column",	&outc,	 	&arg_setint, 	"an error occured", "sets the number of column in the output (default 6)");

	env.optional_arg('C',"incolumn",	&inC,	 	&arg_setint, 	"an error occured", "overides the default number of column assumed for the input (default )");
	env.optional_arg('i',"input",	&infiles,	&arg_setvectorstr,
									"an error occured", "sets the input file (default stdin). If multiple input files given, the default behavior is to merge the files.");
	env.optional_arg('I',"sample",	&sample, 	&arg_setint,	"an error occured", "limit output to sample N.");
	env.optional_arg('o',"output",	&outfile,	&arg_setstr, 	"an error occured", "sets the output file (default stdout)");
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occured while displaying the help message", "prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occured while displaying the version message", "prints the program version");
	env.flag(	'P',"pro",	&args.pro, 	&flag_set, 	"an error occured", "input is in the pro format");
	env.flag(	'b',"binary",	&binary, 	&flag_set, 	"an error occured", "output in a binary format");
	env.flag(	'B',"bernhard",	&bernhard, 	&flag_set, 	"an error occured", "print in mlRho compatibility mode.");
	env.flag(	'H',"noheader",	&noheader, 	&flag_set, 	"an error occured", "don't print a header.");
	

	if (parsargs(argc, argv, env)!=0) exit(0);

	if (bernhard) { 
		outc=5; 
		noheader=true; 
		qdel='\t'; 
		if(sample==-1) sample=1;
	} 

	if ( outc!=6 && outc!=5 && outc!=7 ) {std::cerr << "columns must be 5, 6 or 7 (e.g. -c 5).\n"; exit(0);}

	//Setting up the input/output

	profile out;			//the output profile
	std::vector <profile *> in;	//the input profile(s)

	count_t samples=0;

	if (infiles.size()!=0){
		for (size_t x=0; x<infiles.size(); x++){
			in.push_back(new profile);
			in[x]->open(infiles[x].c_str(), std::fstream::in);
			if (!in[x]->is_open() ) {printUsage(env); exit(0);}
			if (inC!=0) in[x]->setcolumns(inC);
			samples+=in[x]->size();
		}
	} else {
		in.push_back(new profile);
		in[0]->open(std::fstream::in);
		if (!in[0]->is_open() ) {printUsage(env); exit(0);}
		if (inC!=0) in[0]->setcolumns(inC);
		samples+=in[0]->size();
		
	}

	if (binary) {
		if (outfile.size()!=0) {
			out.open(outfile.c_str(), std::fstream::out | std::fstream::binary);
			if (!out.is_open() ) {printUsage(env); exit(0);} 
		}
		else out.open(std::fstream::out | std::fstream::binary);
	} else {
		if (outfile.size()!=0) {
			out.open(outfile.c_str(), std::fstream::out);
			if (!out.is_open() ) {printUsage(env); exit(0);} 
		}
		else out.open(std::fstream::out);
	}

	out.setsamples(samples);
	out.setcolumns(outc);
	out.set_delim_column(cdel);
	out.set_delim_quartet(qdel);

	count_t z=0;

	//change to iterator...

	size_t thissample=0;

	if (sample==-1){
		for (size_t x=0; x<in.size(); x++){		//copies the sample names from the in file(s) to out.
			for (size_t y=0;y<in[x]->size(); ++y){
				out.setsample_name(z, in[x]->getsample_name(y) );
				z++;
			}
		}
	}
	else {
		for (size_t x=0; x<in.size(); x++){		//copies the sample names from the in file(s) to out.
			for (size_t y=0;y<in[x]->size(); ++y){
			++thissample;
			if (thissample==sample)
			out.setsample_name(0, in[x]->getsample_name(y) );
			}
		}
	}

	if (!noheader) out.writeheader();		//writes a header iff noheader is not set.
	else out.noheader();

	out.set_id0(-1);
	out.set_id1(-1);	//set the ids of out to -1 for syncronization purposes.
	
	for (size_t x=0; x<in.size(); x++){ in[x]->read(); }

	bool go=true;
	bool read_site=false;
	bool read_scaffold=false;

	Locus site=out.get_locus();

	while(go){
		std::vector <quartet_t>::iterator it=site.begin();
		std::vector <quartet_t>::iterator end=site.end();
		go=false;
		read_site=false;
		read_scaffold=false;
		for (int x=0; x<in.size(); x++){
			std::vector <quartet_t>::iterator it_in=in[x]->begin();
			std::vector <quartet_t>::iterator end_in=in[x]->end();	
			if (in[x]->get_id1()==site.get_id1() && out.encodeid0(in[x]->decodeid0(in[x]->get_id0()))==site.get_id0() ){
				site.set_extraid(in[x]->get_extraid(0), 0);
				while (it_in!=end_in){
					*it=*it_in;
					it++;
					it_in++;
				}
				read_site=true;
				read_scaffold=true;
				if (in[x]->read()==EOF ) in[x]->close();
			} else {
				if (out.encodeid0(in[x]->decodeid0(in[x]->get_id0()))==site.get_id0() ) read_scaffold=true;
				while (it_in!=end_in){ 
					*it=0;
					it++;
					it_in++;
				}
			} 
			if (in[x]->is_open() ) go=true;
		}
		if (read_site){
			if (bernhard){
				Locus L=site;
				L.resize(1);
				L.set_quartet(site.get_quartet(sample), 0);
				out.write(L);
			} else out.write(site);
		} 
		if (!read_scaffold){
			site.set_id0(site.get_id0()+1);
			site.set_id1(0);
		} else {
			site.set_id1(site.get_id1()+1);
		}
	};

	out.close();
	for (size_t x=0; x<in.size(); ++x) {in[x]->close(); delete in[x];}
	env.close();
	return 0;
}


