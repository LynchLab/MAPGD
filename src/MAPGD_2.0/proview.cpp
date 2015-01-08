/***** proview.c **********************************
 * This is sam2pro almost verbatem.
 * Description: Convert sam output to profiles.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jul 21 22:46:11 2010
 **************************************************/

#include "proview.h"

int proview(int argc, char *argv[])
{
	/* commands that can be set from the command line */

	/* default values for arguments */
	Args args;
	args.delim='\t';
	args.min=4;
	args.c=5;
	args.pvalue=0.001;
	args.notrim=0;

	std::string infile="";
	std::string outfile="";

	std::ifstream inFile;				// the file to read data from (if not stdin).
	std::ofstream outFile;				// the file to write data to (if not stdout).

	env_t env;
	env.setname("mapgd proview");
	env.setver("0.6a");
	env.setauthor("Bernhard Haubold");
	env.setdescription("prints data in the '.pro' file quartet format");
	
	env.optional_arg('m',"minimum",	 &args.min, 	&arg_setint, 	"an error occured", "prints a line iff at least one line has coverage greater than the minimum coverage (defauld 4)");
	env.optional_arg('d',"delimiter",&args.delim,	&arg_setchar, 	"an error occured", "sets the column delimiter (default tab)");
	env.optional_arg('c',"coulmn",	&args.c, 	&arg_setint, 	"an error occured", "sets the number of column in the output (default 6)");
	env.optional_arg('i',"input",	&infile,	&arg_setstr, 	"an error occured", "sets the input file (default stdin)");
	env.optional_arg('o',"output",	&outfile,	&arg_setstr, 	"an error occured", "sets the output file (default stdout)");
	env.optional_arg('t',"trim",	&args.pvalue,	&arg_setdouble,	"an error occured", "skip printing lines where an allele occurs primarly in one direction, \n\t\tgive that the p-value < the number provided");
	env.flag(	'n',"notrim",	&args.notrim,	&flag_set,	"an error occured", "disable trimming");
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occured while displaying the help message", "prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occured while displaying the version message", "prints the program version");

	parsargs(argc, argv, env) ;

	if ( args.c!=6 && args.c!=5 ) {std::cerr << "columns must be 5 and 6 (e.g. -c 5 or -c 6).\n"; exit(0);}

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

	std::ostream *out;
	std::istream *in;

	out=&std::cout;
	in=&std::cin;

	if (infile.size()!=0) {
		inFile.open(infile, std::ifstream::in);
		if (!inFile) {printUsage(env);} 
		else in=&inFile;
	}
	if (outfile.size()!=0) {
		outFile.open(outfile, std::ofstream::out);
		if (!outFile) {printUsage(env);} 
		else out=&outFile;
	}

	/* Open the input file. */

	/* ************************************************************************************************************ */

	runAnalysis(dic, args, in, out);
	if(outfile.size()!=0) outFile.close();
	if(infile.size()!=0) inFile.close();
	env.close();
	exit(0);
}

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



void runAnalysis(int *dic, Args args, std::istream *in, std::ostream *out)
{
	int s, fieldc=0, linec=0;
	char consensus; 
	std::string name, column, number, line;
	std::vector <std::string> field;

	getline(*in, line);
	field=split(line, args.delim);
	unsigned int samples=(field.size()-3)/3;
	int ** count=new int* [samples];

	for(int i = 0 ; i < samples; ++i)
		count[i] = new int[9];
	
	/*If multiple files are given as arguments, the default behavoir should be to merge the files, an option should be given to append the files, 
	* but in the event that duplicates of a site exist the program should throw an error or warning.
	*/

	while( std::getline(*in, line)!=NULL){
		field=split(line, args.delim);
		++linec;
		fieldc = field.size();
		if(fieldc < 5){
			printf("WARNING[mapgd proview]: Skipping line %d with only %d fields.\n", linec, fieldc);
			continue;
		}
		for(int i = 0 ; i < samples; ++i)
    			memset(count[i], 0, sizeof(int)*9 );

		if(name.compare(field[0]) != 0){
			name = field[0];
			if(args.c == 5) (*out) << ">" << name << std::endl;
		}

		consensus = field[2][0];
		for (int i = 0; i < samples; ++i){
			column = field[4+i*3];
			scanCol(column, dic, count[i], consensus);
		}

		bool run=false;
		
		for (int x=0 ; x<samples ; ++x){
			s=0;
			for(int y=0; y<8; ++y) s += count[x][y];
			if (s>=args.min) run=true;
		}
		if (!args.notrim){
			for (int x=0 ; x<samples ; ++x)
				for (int y=0 ; y<4 ; y++) 
					if (biased(count[x][y], count[x][y+4], args.pvalue) ) run=false;
		}
		if(run){
			if(args.c == 6) (*out) << name << "\t";
			(*out) << field[1];
			for(int x=0 ; x<samples ; x++){ 
				for (int y=0 ; y<4 ; y++){ 
					(*out) << "\t" << count[x][y]+count[x][y+4];
			} }
			(*out) << std::endl;
		}
	}
}

void scanCol(std::string column, int *dic, int *count, char consensus)
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
