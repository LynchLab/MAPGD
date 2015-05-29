/***** proview.cpp ********************************
 * This is sam2pro almost verbatem.
 * Description: Convert sam output to profiles.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Sat Apr 11 09:44:00 2014
 **************************************************/

#include "proview.h"

int proview(int argc, char *argv[])
{
	/* commands that can be set from the command line */

	/* default values for arguments */
	Args args;
	args.cdel='\t';
	args.qdel='/';
	args.pro=false;
	args.min=4;
	args.c=5;
	args.pvalue=0.001;
	args.notrim=false;
	args.noheader=false;

	std::string infile="";
	std::string outfile="";

	std::ifstream inFile;				// the file to read data from (if not stdin).

	env_t env;
	env.setname("mapgd proview");
	env.setver("0.6a");
	env.setauthor("Bernhard Haubold");
	env.setdescription("prints data in the '.pro' file quartet format");
	
	env.optional_arg('m',"minimum",	 &args.min, 	&arg_setint, 	"an error occured", "prints a line iff at least one line has coverage greater than the minimum coverage (defauld 4)");
	env.optional_arg('q',"qdel",	&args.qdel,	&arg_setchar, 	"an error occured", "sets the quartet delimiter (default tab)");
	env.optional_arg('d',"cdel",	&args.cdel,	&arg_setchar, 	"an error occured", "sets the column delimiter (default tab)");
	env.optional_arg('c',"coulmn",	&args.c, 	&arg_setint, 	"an error occured", "sets the number of column in the output (default 6)");
	env.optional_arg('i',"input",	&infile,	&arg_setstr, 	"an error occured", "sets the input file (default stdin)");
	env.optional_arg('o',"output",	&outfile,	&arg_setstr, 	"an error occured", "sets the output file (default stdout)");
	env.optional_arg('t',"trim",	&args.pvalue,	&arg_setfloat_t,	"an error occured", "skip printing lines where an allele occurs primarly in one direction, \n\t\tgive that the p-value < the number provided");
	env.flag(	'n',"notrim",	&args.notrim,	&flag_set,	"an error occured", "disable trimming");
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occured while displaying the help message", "prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occured while displaying the version message", "prints the program version");
	env.flag(	'P',"pro",	&args.pro, 	&flag_set, 	"an error occured", "input is a the pro format");
	env.flag(	'N',"noheader",	&args.noheader, &flag_set, 	"an error occured", "don't print a header.");
	parsargs(argc, argv, env) ;

	if ( args.c!=6 && args.c!=5 && args.c!=7) {std::cerr << "columns must be 5, 6 or 7 (e.g. -c 5).\n"; exit(0);}

	int dic[256]; //A dictionary for translating text characters to ints, currently set to map forward reads and reverse reads to different ints.

	//memset(dic, 8, sizeof(int)*256);
	for (int x=0; x<256; ++x)dic[x]=8;

	dic['A'] = 0;
	dic['C'] = 1;
	dic['G'] = 2;
	dic['T'] = 3;
	dic['a'] = 4;
	dic['c'] = 5;
	dic['g'] = 6;
	dic['t'] = 7;

	//Setting up the input/output

	profile pro;
	profile pro_in;
	std::istream *in;
	in=&std::cin;

	if (not (args.pro) ){	
		if (infile.size()!=0) {
			inFile.open(infile, std::ifstream::in);
			if (!inFile) {printUsage(env);} 
			else in=&inFile;
		}
	}
	else{
		if (infile.size()!=0) pro_in.open(infile.c_str(), 'r');	
		else pro_in.open('r');
	};
	if (outfile.size()!=0) {
		pro.open(outfile.c_str(), 'w');
		if (!pro.is_open() ) {printUsage(env); exit(0);} 
	}
	else pro.open('w');

	/* Open the input file. */

	/* ************************************************************************************************************ */

	if (not (args.pro) ) runAnalysis(dic, args, in, pro);
	else {
		pro_in.readheader();
		pro.copyheader(pro_in);
		pro.setcolumns(args.c);
		pro.set_delim_column(args.cdel);
		pro.set_delim_quartet(args.qdel);
		if (not (args.noheader) ) pro.writeheader();

		while(pro_in.read()!=EOF ){
			pro.copy(pro_in);
			pro.write();
		}
	};
	pro.close();
	if(infile.size()!=0){
		if (args.pro) inFile.close();
		else pro_in.close();
	};
	env.close();
	exit(0);
}

/*THIS IS NOT WORKING YET!! FIX IT SOMEDAY!!*/
int biased(int F, int R, double sig)
{
	int m, M;
	if (F<R){ m=F; M=R;}
	else m=R; M=F;
	if (F+R>40){
		if ( m>=( (F+M)/2.0-0.5*(1.0+erf(sig/SQRT2 ) ) ) )
			return false;
		else return true;
	} else {
		return false;
		if (m>=( (F+M)/2.0-3.719*sqrt((F+M)/4.0) ) )
			return false;
		else return true;
	}
}



void runAnalysis(int *dic, Args args, std::istream *in, profile &pro)
{

	int s, fieldc=0, linec=0;
	char consensus;

	std::string name, column, number, line;
	std::vector <std::string> field;


	getline(*in, line);
	field=split(line, args.cdel);

	count_t samples=(field.size()-3)/3;
	count_t ** count=new count_t* [samples];

	for(count_t i = 0 ; i < samples; ++i)
		count[i] = new count_t[9];

	site_t site(samples);
	pro.setsamples(samples);
	pro.setcolumns(args.c);
	pro.set_delim_column(args.cdel);
	pro.set_delim_quartet(args.qdel);
	if (not (args.noheader) ) pro.writeheader();

	/* If multiple files are given as arguments, the default behavoir should be to merge the files,
	*  an option should be given to append the files, but in the event that duplicat ids exist in 
	*  files being appended the program should throw an error or warning.
	*/

	while( std::getline(*in, line)!=NULL){
		field=split(line, args.cdel);
		++linec;
		fieldc = field.size();
		if(fieldc < 5){
			printf("WARNING[mapgd proview]: Skipping line %d with only %d fields.\n", linec, fieldc);
			continue;
		}
		for(count_t i = 0 ; i < samples; ++i)
    			memset(count[i], 0, sizeof(count_t)*9 );

		site.id0 = field[0];
		consensus = field[2][0];
		if (site.extra_ids.size()==0) site.extra_ids.push_back(std::string(1, consensus));
		else site.extra_ids[0] = consensus;

		for (count_t i = 0; i < samples; ++i){
			column = field[4+i*3];
			scanCol(column, dic, count[i], consensus);
		}

		bool run=false;
		
		for (count_t x=0 ; x<samples ; ++x){
			s=0;
			for(count_t y=0; y<8; ++y) s += count[x][y];
			if (s>=args.min) run=true;
		}
		if (!args.notrim){
			for (count_t x=0 ; x<samples ; ++x)
				for (count_t y=0 ; y<4 ; y++) 
					if (biased(count[x][y], count[x][y+4], args.pvalue) ) run=false;
		}
		if(run){
			site.id1 = field[1];
			for(count_t x=0 ; x<samples ; x++){ 
				for (count_t y=0 ; y<4 ; y++){ 
					site.sample[x].base[y]=count[x][y]+count[x][y+4];
				} 	
			}
			pro.write(site);
		}
	}
}

void scanCol(std::string column, int *dic, count_t *count, char consensus)
{
	int j;
	char c;
	std::string number;

	for(int i=0 ; i < column.size() ; ++i){
		c = column[i];
		if(c == '$') continue;
    		else if(c == '^') ++i;
		else if(c == '+' || c == '-'){
			number = "";
			while(isdigit(column[++i])) number.push_back(column[i]);
 			for(j=0;j<atoi(number.c_str() )-1;j++) i++;
    		} else if(c == ',' ) {
			count[(int)dic[(int)consensus]]++;
		} else if (c == '.' ) {
			count[(int)dic[(int)consensus]+4]++;
		} else count[(int)dic[(int)c]]++;
	}
}
