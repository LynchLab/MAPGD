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

	count_t outc=6;
	count_t inC=0;

	args.pro=false;
	args.min=4;
	args.pvalue=0.001;

	std::vector <std::string> infiles;
	std::string outfile="";

	env_t env;
	env.setname("mapgd proview");
	env.setver("0.6a");
	env.setauthor("Bernhard Haubold");
	env.setdescription("prints data in the '.pro' file quartet format");
	
	env.optional_arg('m',"minimum",	 &args.min, 	&arg_setint, 	"an error occured", "prints a line iff at least one line has coverage greater than the minimum coverage (defauld 4)");
	env.optional_arg('q',"qdel",	&qdel,		&arg_setchar, 	"an error occured", "sets the quartet delimiter (default tab)");
	env.optional_arg('d',"cdel",	&cdel,		&arg_setchar, 	"an error occured", "sets the column delimiter (default tab)");
	env.optional_arg('c',"column",	&outc,	 	&arg_setint, 	"an error occured", "sets the number of column in the output (default 6)");

	env.optional_arg('C',"incolumn",	&inC,	 	&arg_setint, 	"an error occured", "overides the default number of column assumed for the input (default )");
	env.optional_arg('i',"input",	&infiles,	&arg_setvectorstr,
									"an error occured", "sets the input file (default stdin). If multiple input files given, the default behavior is to merge the files.");
	env.optional_arg('o',"output",	&outfile,	&arg_setstr, 	"an error occured", "sets the output file (default stdout)");
//	env.optional_arg('t',"trim",	&args.pvalue,	&arg_setfloat_t,	"an error occured", "skip printing lines where an allele occurs primarly in one direction, \n\t\tgive that the p-value < the number provided");
//	env.flag(	'n',"notrim",	&args.notrim,	&flag_set,	"an error occured", "disable trimming");
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occured while displaying the help message", "prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occured while displaying the version message", "prints the program version");
	env.flag(	'P',"pro",	&args.pro, 	&flag_set, 	"an error occured", "input is in the pro format");
	env.flag(	'b',"binary",	&binary, 	&flag_set, 	"an error occured", "output in a binary format");
	env.flag(	'H',"noheader",	&noheader, 	&flag_set, 	"an error occured", "don't print a header.");

	if (parsargs(argc, argv, env)!=0) exit(0);

	if ( outc!=6 && outc!=5 && outc!=7 ) {std::cerr << "columns must be 5, 6 or 7 (e.g. -c 5).\n"; exit(0);}

	//Setting up the input/output

	profile out;			//the output profile
	std::vector <profile *> in;	//the input profile(s)

	count_t samples=0;

	if (infiles.size()!=0){
		for (int x=0; x<infiles.size(); x++){
			in.push_back(new profile);
			in[x]->open(infiles[x].c_str(), "r");
			if (!in[x]->is_open() ) {printUsage(env); exit(0);}
			if (inC!=0) in[x]->setcolumns(inC);
			//in[x]->set_delim_column(cdel);
			//in[x]->set_delim_quartet(qdel);
			samples+=in[x]->size();
		}
	} else {
		in.push_back(new profile);
		in[0]->open("r");
		if (!in[0]->is_open() ) {printUsage(env); exit(0);}
		if (inC!=0) in[0]->setcolumns(inC);
		//in[0]->set_delim_column(cdel);
		//in[0]->set_delim_quartet(qdel);
		samples+=in[0]->size();
	}

	if (binary) {
		if (outfile.size()!=0) {
			out.open(outfile.c_str(), "wb");
			if (!out.is_open() ) {printUsage(env); exit(0);} 
		}
		else out.open("wb");
	} else {
		if (outfile.size()!=0) {
			out.open(outfile.c_str(), "w");
			if (!out.is_open() ) {printUsage(env); exit(0);} 
		}
		else out.open("w");
	}

	out.setcolumns(outc);
	out.setsamples(samples);
	out.set_delim_column(cdel);
	out.set_delim_quartet(qdel);

	count_t z=0;

	//change to iterator...

	for (int x=0; x<in.size(); x++){		//copies the sample names from the in file(s) to out.
		for (int y=0;y<in[x]->size(); ++y){
			out.setsample_name(z, in[x]->getsample_name(y) );
			z++;
		}
	} 

	if (!noheader) out.writeheader();		//writes a header iff noheader is not set.

	out.setid0(0);
	out.setid1(0);					//set the ids of out to 0 for syncronization purposes.

	for (int x=0; x<in.size(); x++) in[x]->read();	//read the first line of each file.
	//in[x]->decodeid(0)

		

	bool go=true;
	while(go){
		go=false;
		std::vector <quartet_t>::iterator it=out.begin();
		std::vector <quartet_t>::iterator end=out.end();
		bool read=false;
		bool wrote=false;
		for (int x=0; x<in.size(); x++){
			if (in[x]->is_open() ){
				if(out.encodeid0(in[x]->decodeid0(in[x]->getid0()))==out.getid0() ){
					if (in[x]->getid1()==out.getid1() ){
						if (in[x]->read()==EOF ){ 
							in[x]->close();
						}
					}
					read=true;
				};
				go=true;
			}

		}
			
		if (!read){
			out.setid0(out.getid0()+1);
			out.setid1(0);
		} else {
			out.setid1(out.getid1()+1);
		}
		for (int x=0; x<in.size(); x++){
			std::vector <quartet_t>::iterator it_in=in[x]->begin();
			std::vector <quartet_t>::iterator end_in=in[x]->end();	
			if (in[x]->getid1()==out.getid1() && out.encodeid0(in[x]->decodeid0(in[x]->getid0()))==out.getid0() ){
				out.setextraid(in[x]->getextraid(0), 0);
				//change to memcopy;
				wrote=true;
				while (it_in!=end_in){ 
					*it=*it_in;
					it++;
					it_in++;
				}
			}
			else{
				//change to memset;
				while (it_in!=end_in){ 
					*it=0;
					it++;
					it_in++;
				}
			}
		}
		if (wrote) out.write();
	};

	out.close();
	for (int x=0; x<in.size(); ++x) {in[x]->close(); delete in[x];}

	env.close();
	exit(0);
}


