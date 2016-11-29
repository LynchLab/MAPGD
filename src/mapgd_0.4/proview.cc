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

//class Bcf2pro : public Locus{};
//class Bcf2pro_file : public Indexed_file <Bcf2pro>{};


using namespace std;

int proview(int argc, char *argv[])
{
	/* commands that can be set from the command line */

	/* default values for arguments */
	Args args;

	bool out_binary=false;
	bool in_pro=false;
	bool bernhard=false;
	bool noheader=false;
	bool dontprint=false;

	args.pro=false;
	args.min=4;
	args.pvalue=0.001;

	std::vector <std::string> infiles;
	std::string namefile="";
	std::string outfile="";
	std::string headerfile="";

	Region region(0,ID1_MAX);

	Environment env;
	env.set_name("mapgd proview");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman and Bernhard Haubold");
	env.set_description("prints data in the '.pro' file quartet format");
	env.required_arg('H',"header",	headerfile,	"You must specify an index file (-H)", "sets the index file (required to use mpileup)");
	env.optional_arg('m',"minimum",	args.min, 	"please provide a integer number", "prints a line iff at least one line has coverage greater than the minimum coverage (defauld 4)");
	env.optional_arg('R',"region",	region, 	"please provide a valid region (name:start-stop)", "a string which specifies the region to be viewed");
	env.optional_arg('n',"names",	namefile, 	"No name file specified", "a tab delimited file with sample name 'tab' file name pairs");
	env.optional_arg('o',"output",	outfile,	"No output file specified", "sets the output file (default stdout)");
	env.positional_arg('i',"input",	infiles,	"No input file specified", "the mpileup files to be used");
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message", "prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occurred while displaying the version message", "prints the program version");
	env.flag(	'b',"binary",	&out_binary, 	&flag_set, 	"an error occurred", "output in a binary format");
	env.flag(	's',"skip",	&dontprint,	&flag_set, 	"an error occurred", "skip empty lines");
	env.flag(	'r',"mlrho",	&bernhard, 	&flag_set, 	"an error occurred", "output in mlrho format");
	env.flag(	'N',"noheader",	&noheader,	&flag_set, 	"an error occurred", "don't print those silly '@' lines. Overrides the b option.");
	env.flag(	'p',"pro",	&in_pro, 	&flag_set, 	"an error occurred", "input is in pro format");

	if (parsargs(argc, argv, env)!=0) exit(0);
	
	Indexed_file <Locus> out_file;		//the output profile

	std::vector <Bcf2pro_file *> in_files;	//the input profile(s)

	std::vector <Locus> in_locus;
	Locus out_locus;

	File_index index;

        if (headerfile.size()!=0) {
                std::fstream header;
                header.open(headerfile.c_str(),  ios::in);
                index.from_sam_header( header );
                out_file.set_index(index);
        }

	if (index.get_sizes().size()==0) {
		std::cerr << __FILE__ << ":" << __LINE__ << ". Error: no scaffolds in index file. Exiting.\n";
		exit(0);
	}

	std::vector <std::string> sample_names;

	size_t sample_numbers=0;

	if (namefile.size()!=0) {
		Flat_file <Sample_name> in_names;
		Sample_name name_file;
		in_names.open(namefile.c_str(), ios::in);
		name_file=in_names.read_header();
		while (in_names.read(name_file).table_is_open() ){
			in_files.push_back(new Bcf2pro_file(in_pro) );
			in_files.back()->open_no_extention(name_file.mpileup_name.c_str(), ios::in);
			in_locus.push_back( in_files.back()->read_header() );
			if (name_file.sample_names.size()!=in_locus.back().get_sample_names().size() ){
				std::cerr << __FILE__ << ":" << __LINE__ << ". Error: name file does not name the correct number of samples. Exiting.\n";
				exit(0);
			}
			for (size_t y=0; y<name_file.sample_names.size(); ++y){
				sample_names.push_back(name_file.sample_names[y]);
			}
			sample_numbers+=in_locus.back().get_sample_names().size();
       		 	in_files.back()->set_index(index);
			if (in_files.back()->read(in_locus.back()).eof() ) in_files.back()->close();
		}
		in_names.close();
	} else	if (infiles.size()!=0){	
		for (size_t x=0; x<infiles.size(); ++x) {
			in_files.push_back(new Bcf2pro_file(in_pro) );
			if (infiles[x].size()!=0) in_files.back()->open_no_extention(infiles[x].c_str(), ios::in);
			else in_files.back()->open(ios::in);
			in_locus.push_back( in_files.back()->read_header() );
			for (size_t y=0; y<in_locus.back().get_sample_names().size(); ++y){
				std::stringstream s;
				s << split_last(infiles[x], '/').back() << ":" << y+1;
				sample_names.push_back( s.str().c_str() );
			}
			sample_numbers+=in_locus.back().get_sample_names().size();
              	 	in_files[x]->set_index(index);
			if (in_files[x]->read(in_locus[x]).eof() ) in_files[x]->close();
		}
	} else 	{
		in_files.push_back(new Bcf2pro_file(in_pro) );
		in_files.back()->open(ios::in);
		in_locus.push_back( in_files.back()->read_header() );
		for (size_t y=0; y<in_locus.back().get_sample_names().size(); ++y){
			sample_names.push_back(in_locus.back().get_sample_names()[y]);
		}
                in_files[0]->set_index(index);
		sample_numbers+=in_locus.back().get_sample_names().size();
		if (in_files[0]->read(in_locus[0]).eof() ) in_files[0]->close();
	}

	if (sample_numbers==0) {
		std::cerr << __FILE__ << ":" << __LINE__ << ". Error: no mpileup files opened. Exiting.\n";
		std::cerr << "You can generate mpileup files with the command samtools mpileup.\n";
		exit(0);
	}

	if (outfile.size()==0) out_file.open( out_binary ? ios::out | ios::binary : ios::out );
	else out_file.open(outfile.c_str(), out_binary ? ios::out | ios::binary : ios::out );
	out_locus.set_sample_names(sample_names);
	out_file.set_index(index);

	if(!noheader) {
		out_file.write_header(out_locus);
	} else {
	//	out_file.write_header(out_locus);
	}

	out_locus.set_abs_pos(1);
	bool print_all=!dontprint, go=true, read_site=false;

	region.set(index);
	id1_t out_abs_pos;

	while(go){
		out_locus.ref.base=4;

		std::vector <quartet_t>::iterator it=out_locus.begin();
		std::vector <quartet_t>::iterator end=out_locus.end();

		read_site=false;
		
		for (size_t x=0; x<in_files.size(); x++){

			std::vector <quartet_t>::iterator it_in=in_locus[x].begin();
			std::vector <quartet_t>::iterator end_in=in_locus[x].end();	

			out_locus.ref.base=4;
			out_file.get_pos(out_locus);
			in_files[0]->get_pos(in_locus[0]);
			if ( in_files[x]->is_open()  && in_files[x]->get_pos(in_locus[x])==out_file.get_pos(out_locus) ){
				read_site=true;
				if(in_locus[x].ref.base!=4) out_locus.ref=in_locus[x].ref.base;
				while (it_in!=end_in){
					*it=*it_in;
					it++;
					it_in++;
				}
				if (in_files[x]->read(in_locus[x]).eof() ){
					 in_files[x]->close();
				}
			} else {
				while (it_in!=end_in){ 
					*it=0;
					it++;
					it_in++;
				}
			}
		}
		out_abs_pos=out_locus.get_abs_pos();
		if ( read_site || print_all ){
			if ( out_abs_pos > region.abs_start && out_abs_pos < region.abs_stop ){
				out_locus.unmaskall();
				out_file.write(out_locus);
			}
		}
		out_locus.set_abs_pos(out_abs_pos+1);
		if (out_abs_pos<=region.abs_stop ) go=true;
		else go=false;
	
	};
	if(!noheader) out_file.close();
	for (size_t x=0; x<in_files.size(); ++x) {
		in_files[x]->close(); 
		delete in_files[x];
	}
	env.close();
	return 0;
}


