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
	bool out_cov=false;
	
	int offset=2;
	int columns=3;
    count_t GC=0, N=0, K=0;

	args.pro=false;
	args.min=4;
	args.pvalue=0.001;

	std::vector <std::string> infiles;
	std::vector <std::string> name_list;
	std::string namefile="";
	std::string outfile="";
	std::string headerfile="";

	Region region(0,ID1_MAX);
	std::vector <size_t> ind;

	Environment env;
	env.set_name("mapgd proview");
	env.set_version(VERSION);
	env.set_author("Matthew Ackerman and Bernhard Haubold");
	env.set_description("prints data in the '.pro' file quartet format");
	env.optional_arg('I',"individuals", ind, 	"please provide a list of integers", "the individuals to be present in the output profile. A comma separated list containing no spaces, the python slice notation (i.e. 1:4,7 for 1,2,3,4,7 or 1:5:2,7 for 1,3,5,7 ... can be used to specify this list (default ALL).");
	env.optional_arg('H',"header",	headerfile,	"You must specify an index file (-H)", "sets the index file (required to use mpileup)");
	env.optional_arg('m',"minimum",	args.min, 	"please provide a integer number", "prints a line iff at least one line has coverage greater than the minimum coverage (defauld 4)");
	env.optional_arg('R',"region",	region, 	"please provide a valid region (name:start-stop)", "a string which specifies the region to be viewed");
	env.optional_arg('f',"offset",	offset, 	"please provide a valid integer", "offset untill the first sample column in mpileup");
	env.optional_arg('c',"columns",	columns, 	"please provide a valid integer", "number of columns per sample");
	env.optional_arg('n',"names",	namefile, 	"No name file specified", "a tab delimited file with sample name 'tab' file name pairs");
	env.optional_arg('l',"name_list",	name_list, 	"please provide a comma separated list", "a comma delimited list used as names for the samples.");
	env.optional_arg('o',"output",	outfile,	"No output file specified", "sets the output file (default stdout)");
	env.positional_arg('i',"input",	infiles,	"No input file specified", "the mpileup files to be used");
	env.flag(	'h',"help", 	&env, 		&flag_help, 	"an error occurred while displaying the help message", "prints this message");
	env.flag(	'v',"version", 	&env, 		&flag_version, 	"an error occurred while displaying the version message", "prints the program version");
	env.flag(	'b',"binary",	&out_binary, 	&flag_set, 	"an error occurred", "output in a binary format");
	env.flag(	's',"skip",	&dontprint,	&flag_set, 	"an error occurred", "skip empty lines");
	env.flag(	'r',"mlrho",	&bernhard, 	&flag_set, 	"an error occurred", "output in mlrho format");
	env.flag(	'N',"noheader",	&noheader,	&flag_set, 	"an error occurred", "don't print those silly '@' lines. Overrides the b option.");
	env.flag(	'p',"pro",	&in_pro, 	&flag_set, 	"an error occurred", "input is in pro format");
	env.flag(	'd',"depth",	&out_cov, 	&flag_set, 	"an error occurred", "output coverage instead of profile.");

	if (parsargs(argc, argv, env)!=0)
	{
	//	if (!in_pro) 
		exit(0);
	}
	if( headerfile=="" and !in_pro)
	{
	    fprintf(stderr, gettext("mapgd:%s:%d: You must either specify an index file (-H) or read in a profile (-p).\n"), __FILE__, __LINE__);
	}
	
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

	if (!in_pro && index.get_sizes().size()==0) {
	    fprintf(stderr, gettext("mapgd:%s:%d: no scaffolds in index file.\n"), __FILE__, __LINE__);
		exit(0);
	}

	std::vector <std::string> sample_names;

	size_t sample_numbers=0;

	if (namefile.size()!=0) {
        if (name_list.size()!=0) {
	        fprintf(stderr, gettext("mapgd:%s:%d: Warning: a name list (-l) and a name file (-n) have both been specified. Ignoring name list.\n"), __FILE__, __LINE__);
        }
		Flat_file <Sample_name> in_names;
		Sample_name name_file;
		in_names.open(namefile.c_str(), ios::in);
		name_file=in_names.read_header();
		while (in_names.read(name_file).table_is_open() ){
			in_files.push_back(new Bcf2pro_file(in_pro) );
            if (name_file.mpileup_name.c_str()[0]=='-')
            {
    			in_files.back()->open(ios::in);
            } else {
    			in_files.back()->open_no_extention(name_file.mpileup_name.c_str(), ios::in);
            }
            in_files.back()->set_mpileup(offset, columns);
			in_locus.push_back( in_files.back()->read_header() );
			if (name_file.sample_names.size()!=in_locus.back().get_sample_names().size() ){
	            fprintf(stderr, gettext("mapgd:%s:%d: name file %s does not name the correct number of samples.\n"), __FILE__, __LINE__, name_file.mpileup_name.c_str() );
				exit(0);
			}
			for (size_t y=0; y<name_file.sample_names.size(); ++y){
				sample_names.push_back(name_file.sample_names[y]);
			}
			sample_numbers+=in_locus.back().get_sample_names().size();
       		 	if (!in_pro) in_files.back()->set_index(index);
			if (in_files.back()->read(in_locus.back()).eof() ) in_files.back()->close();
		}
		in_names.close();
	} else {
        if (name_list.size()!=0) {
            sample_names=name_list;
        }
        if (infiles.size()!=0){	
            for (size_t x=0; x<infiles.size(); ++x) {
                in_files.push_back(new Bcf2pro_file(in_pro) );
                if (infiles[x].size()!=0) in_files.back()->open_no_extention(infiles[x].c_str(), ios::in);
                else in_files.back()->open(ios::in);
                in_files.back()->set_mpileup(offset, columns);
                in_locus.push_back( in_files.back()->read_header() );
               
                if (!in_pro) {
                    for (size_t y=0; y<in_locus.back().get_sample_names().size(); ++y){
                        std::stringstream s;
                        s << split_last(infiles[x], '/').back() << ":" << y+1;
                        sample_names.push_back( s.str().c_str() );
                    }
                } else {
                    for (size_t y=0; y<in_locus.back().get_sample_names().size(); ++y){
                        sample_names.push_back(in_locus.back().get_sample_names()[y]);
                    }
                }
                sample_numbers+=in_locus.back().get_sample_names().size();
                if (!in_pro) in_files[x]->set_index(index);
                if (in_files[x]->read(in_locus[x]).eof() ) in_files[x]->close();
            }
        } else 	{
            in_files.push_back(new Bcf2pro_file(in_pro) );
            in_files.back()->open(ios::in);
            in_files.back()->set_mpileup(offset, columns);
            in_locus.push_back( in_files.back()->read_header() );
            for (size_t y=0; y<in_locus.back().get_sample_names().size(); ++y){
                sample_names.push_back(in_locus.back().get_sample_names()[y]);
            }
            if (!in_pro) in_files[0]->set_index(index);
            sample_numbers+=in_locus.back().get_sample_names().size();
            if (in_files[0]->read(in_locus[0]).eof() ) in_files[0]->close();
        }
    }

	if (sample_numbers==0) {
		std::cerr << __FILE__ << ":" << __LINE__ << ". Error: no mpileup files opened. Exiting.\n";
		std::cerr << "You can generate mpileup files with the command samtools mpileup.\n";
		exit(0);
	}

	if (outfile.size()==0) out_file.open( out_binary ? ios::out | ios::binary : ios::out );
	else out_file.open(outfile.c_str(), out_binary ? ios::out | ios::binary : ios::out );
	out_locus.set_sample_names(sample_names);

	index=in_files.back()->get_index();

	out_file.set_index(index);

	out_locus.redact_all();

	if (ind.size()==0)
		for (int x=0; x<sample_names.size(); ++x)
			ind.push_back(x);

	for (int x=0; x<ind.size(); ++x) {
		out_locus.unredact(ind[x]);
	}

	if(!noheader) {
		out_file.write_header(out_locus);
	} else {
	//	out_file.write_header(out_locus);
	}

	out_locus.set_abs_pos(1);
	bool print_all=!dontprint, go=true, all_closed=true, read_site=false;

	region.set(index);
	id1_t out_abs_pos;

	while(go){
		out_locus.ref.base=4;

		std::vector <quartet_t>::iterator it=out_locus.begin();
		std::vector <quartet_t>::iterator end=out_locus.end();

		read_site=false;
        all_closed=true;
		
		for (size_t x=0; x<in_files.size(); x++){

			std::vector <quartet_t>::iterator it_in=in_locus[x].begin();
			std::vector <quartet_t>::iterator end_in=in_locus[x].end();	

			//out_locus.ref.base=4;
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
				if ( in_files[x]->read(in_locus[x]).eof() ){
					 in_files[x]->close();
				}
			} else {
				while (it_in!=end_in){ 
					*it=0;
					it++;
					it_in++;
				}
            }
            if ( in_files[x]->table_is_open() ) {
                all_closed=false;
                if (in_files[x]->get_pos(in_locus[x])<out_file.get_pos(out_locus) ){
                    if (in_files[x]->get_pos(in_locus[x])<out_file.get_pos(out_locus) ){
                        if (!all_closed) {
                            fprintf(stderr, gettext("mapgd:%s:%d: Error syncing input and output file.\n"), __FILE__, __LINE__);
                            exit(0);
                        }
                    }
               }
           }
		}
		out_abs_pos=out_locus.get_abs_pos();
		if ( read_site || print_all ){
			if ( out_abs_pos > region.abs_start && out_abs_pos < region.abs_stop ){
				out_locus.unmaskall();
                if (!out_cov)
                    out_file.write(out_locus);
                else {
                    {
                        std::cout.width(14);
                        std::cout << std::left << index.get_string(index.get_id0(out_locus.get_abs_pos()) );
                        std::cout << '\t';
                        std::cout << index.get_id1(out_locus.get_abs_pos());
                        std::cout << '\t';
                        std::cout << out_locus.getcoverage();
                        std::cout << '\t';
                        std::cout << int(out_locus.ref=='G' or out_locus.ref=='C');
                        std::cout << '\t';
                        std::cout << int(not(out_locus.ref=='N')) << std::endl;
                    }
                }
			}
		}
		out_locus.set_abs_pos(out_abs_pos+1);
		if (out_abs_pos <= region.abs_stop ) go=true;
		else go=false;
        if (!print_all && all_closed) go=false;
	};

    //Print mean coverages summary statistics?
	if(!noheader) {
		out_file.close_table();

        /*Flat_file <Sample_cov> smp_file;
		smp_file.open_from(out_file);
		smp_file.write_header( Sample_gof() );
		for (size_t x=0; x<sample_names.size(); ++x) {
			smp_file.write(Sample_gof(x+1, sample_names[ind[x]], 0 ) );
		}
		smp_file.close();*/
		out_file.close();
	}


	for (size_t x=0; x<in_files.size(); ++x) {
		in_files[x]->close(); 
		delete in_files[x];
	}
	env.close();

    fprintf(stderr, gettext("mapgd:%s:%d: Clean exit. goodbye!\n"), __FILE__, __LINE__);

	return 0;
}


