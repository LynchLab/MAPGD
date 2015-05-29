#include "proFile.h"

/*Design notes: This library is probably poorly thought out. More or less everything works through the 'profile'. However, 
  It may make a whole lot more sense for the 'profile' to be a file like structure that just co-ordinates reading, writing, 
  formating and all that jazz. Most of the work should probably be done by the "site" structure.
 */

/** @brief Read a .pro file in the five column format.
  * @returns 0 iff successful
**/

const std::string profile::names_ ="ACGT*";

site_t::site_t(void){};
site_t::site_t(count_t size){
	sample.assign(size, quartet() );
};

site_t & site_t::operator =(const site_t& arg){
	
        sample=arg.sample; 	
        id0=arg.id0, id1=arg.id1;  
        extra_ids=arg.extra_ids;   
};

int profile::seek(std::streampos pos) {
	/*if (column==5 || in==NULL){
		std::cerr << "seek is currently only implemented on random access files opened in 6 or 7 column mode. Try re-running ei w/o the -s option.\n";
	};*/
	in->seekg(pos);
	read();
	return 0;
}

int profile::setsamples(count_t samples) {
	samples_=samples;
	column_names.clear();
	if (columns_==6 or columns_==7) column_names.push_back("scaffold");
	column_names.push_back("pos");
	if (columns_==7) column_names.push_back("ref");

	for (unsigned int x=0; x<samples_; ++x){
		column_names.push_back("sample_"+std::to_string( (unsigned long long int)(x+1) ) );
		site_.sample.push_back(quartet_t() );
	};
	return NONE;
}

int profile::setcolumns(count_t x) {
	if (column_names.size()==0) setsamples(0);
	if (columns_==5 and x==6){column_names.insert(column_names.begin(), "scaffold");}
	if (columns_==5 and x==7){column_names.insert(column_names.begin(), "scaffold"); column_names.insert(column_names.begin()+2, "ref");}
	if (columns_==6 and x==5){column_names.erase(column_names.begin() );}
	if (columns_==6 and x==7){column_names.insert(column_names.begin()+2, "ref");}
	if (columns_==7 and x==5){column_names.erase(column_names.begin() ); column_names.erase(column_names.begin()+1);}
	if (columns_==7 and x==6){column_names.erase(column_names.begin()+2);}
	columns_=x;
	return NONE;
}

int profile::setcolumn_name(count_t x, std::string str) {
	column_names[x]=str;
	return NONE;
}

int profile::writeheader(void){
	if (out==NULL){
		std::cerr << "attempted writeheader when no outstream was open." << std::endl;
		exit(UNEXPECTED);
	};
	//Headers should contain the following feilds:
	//The version number of proview used to create the file,
	//The number of columns in the file
	//The number of ids in the files
	//Header lines begin with the '@' character (a nod to samtools).
	std::string line;
	std::vector <std::string> column;
	bool notdone_=true;
	while(notdone_){
		*out << "@PR" << delim_column << "VN:" << VER << delim_column <<"CN:" << columns_ << delim_column << "SN:" << samples_ << std::endl;
		*out << "@ID";
		for (unsigned int x=0; x<column_names.size(); ++x)
			*out << delim_column << column_names[x];
		*out << std::endl;
		notdone_=false;
	}
	return NONE;
}

int profile::copyheader(profile const &pro){
	columns_=pro.columns_;
	setsamples(pro.samples_);
	for (unsigned int x=0; x<pro.column_names.size(); ++x){
		column_names[x]=pro.column_names[x];
	};
};

//TODO FIX THIS, IT'S UGLY.

#define H_PARAM	1
#define H_ID	2
#define H_VER	3
#define H_COL	4
#define H_ERR	5
#define H_SAM	6

int hash(std::string str){
	if (str=="@PR") return H_PARAM;
	if (str=="@ID") return H_ID;
	if (str=="VN") return H_VER;
	if (str=="CN") return H_COL;
	if (str=="SN") return H_SAM;
	else return H_ERR;
};

int profile::readheader(void){
	//Headers should contain the following feilds:
	//The version number of proview used to create the file,
	//The number of columns in the file
	//The number of ids in the files
	//Header lines begin with the '@' character (a nod to samtools).
	if (in==NULL){
		std::cerr << "attempted readheader when no instream was open." << std::endl;
		exit(UNEXPECTED);
	};
	std::string line;
	std::vector <std::string> column, args, arg;
	bool notdone_=true, reading_pr=false, reading_id=false;
	while(notdone_){
		//"VN" version number
		//"@IDs" Names of fields.
		if (std::getline(*in, line)!=NULL){
			if (line[0]=='s'){
				std::cerr << "Depricated header format. Please construct pro files with the \'mapgd proview \' command." << std::endl;
				column=split(line, delim_column);
				columns_=7;
				setsamples(column.size()-3);
				for (unsigned int x=0; x<column.size(); ++x) column_names[x]=column[x];
				notdone_=false;
			}
			else if (line[0]=='@'){
				args=split(line, delim_column);
				if (args.size()==0) return BADHEADER;
				for(std::vector<std::string>::iterator argit = args.begin(); argit != args.end(); ++argit){
					arg=split(*argit, ':');
					switch ( hash(arg[0]) ) {
						case H_PARAM:
							reading_pr=true;
							reading_id=false;
							break;
						case H_VER:
							break;
						case H_COL:
							columns_=atoi(arg[1].c_str() );
							break;
						case H_SAM:
							samples_=atoi(arg[1].c_str() );
							setsamples(samples_);
							break;
						case H_ID:
							column_names.clear();
							reading_pr=false;
							reading_id=true;
							notdone_=false;
							break; 
						default :
							if (reading_id) column_names.push_back(*argit);
							else {
								std::cerr << "Warning: unexpected field encountered in header (" << *argit << "). File may not have opened correctly. Try removing this field and running the program again." << std::endl;
								return BADHEADER;
							}
					}
				}
			}
			else{
				std::cerr << "Warning: File lacks header, assuming a five column format produced by the stand-alone program sam2pro. In the future please use mapgd proview to generate pro files." << std::endl;
				columns_=5;
				setsamples(1);
				delim_quartet='\t';
				column_names[0]="ref";
				column_names[1]="bp";
				column_names[2]="sample0";
				return BADHEADER;
			};
		};
	}
	return NONE;
}

int profile::read(){
	read(0);	
};


int profile::read(int arg){
	std::string line;
	//The order of nucleotides read from a quartet file. Replace this with a memcopy of a constant.
	sorted_[0]=0; sorted_[1]=1; sorted_[2]=2; sorted_[3]=3; sorted_[4]=4;

	std::vector <std::string> column, quartet;
	if (in==NULL){
		std::cerr << "attempted read when no instream was open." << std::endl;
		exit(UNEXPECTED);
	};

	if (std::getline(*in, line)!=NULL){
		if (arg==SKIP) return 0;
		//READ A PRO FILE
		if (columns_==5){
			/* check for start of new scaffold*/
			if (line[0]!='>'){
				column=split(line, delim_column);
				if ( (column.size()-1)!=samples_){
					std::cerr << "could not parse line : \"" << line << "\"" << std::endl;
					std::cerr << column.size() << " fields found." << std::endl;
					std::cerr << "delimiter : \'" << delim_column << "\'" << std::endl;
					exit(0);
					return EOF;
				}
				site_.id1=column[0];
				for (unsigned int x=0; x<samples_; ++x){
					quartet=split(column[x+1], delim_quartet);
					site_.sample[x].base[0]=atoi(quartet[0].c_str() );
					site_.sample[x].base[1]=atoi(quartet[1].c_str() );
					site_.sample[x].base[2]=atoi(quartet[2].c_str() );
					site_.sample[x].base[3]=atoi(quartet[3].c_str() );
				};
			} else {
				line.erase(line.begin());
				site_.id0=line;
				read(arg);
			}
		} else if (columns_==6){
			column=split(line, delim_column);
			if ( (column.size()-2)!=samples_){
				std::cerr << "could not parse line : \"" << line << "\"" << std::endl;
				std::cerr << column.size() << " fields found." << std::endl;
				std::cerr << "delimiter : \'" << delim_column << "\'" << std::endl;
				exit(0);
				return EOF;
			};
			site_.id0=column[0];
			site_.id1=column[1];
			for (unsigned int x=0; x<samples_; ++x){
				quartet=split(column[x+2], delim_quartet);
				site_.sample[x].base[0]=atoi(quartet[0].c_str() );
				site_.sample[x].base[1]=atoi(quartet[1].c_str() );
				site_.sample[x].base[2]=atoi(quartet[2].c_str() );
				site_.sample[x].base[3]=atoi(quartet[3].c_str() );
			};
		} else if (columns_==7){
			column=split(line, delim_column);
			if ( (column.size()-3)!=samples_){
				std::cerr << "could not parse line : \"" << line << "\"" << std::endl;
				std::cerr << column.size() << " fields found." << std::endl;
				std::cerr << "delimiter : \'" << delim_column << "\'" << std::endl;
				exit(0);
				return EOF;
			};
			site_.id0=column[0];
			site_.id1=column[1];
			if (site_.extra_ids.size()==0) site_.extra_ids.push_back("");
			site_.extra_ids[0]=column[2];
			for (unsigned int x=0; x<samples_; ++x){
				quartet=split(column[x+3], delim_quartet);
				//std::cout << column[x+3] << ", " << quartet[0] << "/" << quartet[1] << "/" << quartet[2] << "/" << quartet[3] << std::endl;
				site_.sample[x].base[0]=atoi(quartet[0].c_str() );
				site_.sample[x].base[1]=atoi(quartet[1].c_str() );
				site_.sample[x].base[2]=atoi(quartet[2].c_str() );
				site_.sample[x].base[3]=atoi(quartet[3].c_str() );
			};
		}
		return 0;
	}
	return EOF;
};

int profile::copy(const profile &pro){
	site_=pro.site_;
};

int profile::write(void){
	write(site_);
};

int profile::write(const site_t &thissite){
	if (out==NULL){
		std::cerr << "attempted write when no outstream was open." << std::endl;
		exit(UNEXPECTED);
	};
	if (columns_==5){
		if (thissite.id0!=site_.id0){
			*out << '>' << thissite.id0 << std::endl;
		};
		*out << thissite.id1;
		for (unsigned int x=0; x<samples_; ++x){
			if (thissite.sample[x].masked){
				*out << delim_column;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0;
			}
			else {
				*out << delim_column;
				*out << thissite.sample[x].base[0] << delim_quartet;
				*out << thissite.sample[x].base[1] << delim_quartet;
				*out << thissite.sample[x].base[2] << delim_quartet;
				*out << thissite.sample[x].base[3];
			}
		};
		*out << std::endl;
		site_=thissite;
	} else if (columns_==6){
		*out << thissite.id0 << delim_column << thissite.id1;
		for (unsigned int x=0; x<samples_; ++x){
			if (thissite.sample[x].masked){
				*out << delim_column;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0;
			}
			else {
				*out << delim_column;
				*out << thissite.sample[x].base[0] << delim_quartet;
				*out << thissite.sample[x].base[1] << delim_quartet;
				*out << thissite.sample[x].base[2] << delim_quartet;
				*out << thissite.sample[x].base[3];
			}
		};
		*out << std::endl;	
		site_=thissite;
	} else if (columns_==7){
		*out << thissite.id0 << delim_column << thissite.id1 << delim_column << thissite.extra_ids[0];
		for (unsigned int x=0; x<samples_; ++x){
			if (thissite.sample[x].masked){
				*out << delim_column;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0;
			}
			else {
				*out << delim_column;
				*out << thissite.sample[x].base[0] << delim_quartet;
				*out << thissite.sample[x].base[1] << delim_quartet;
				*out << thissite.sample[x].base[2] << delim_quartet;
				*out << thissite.sample[x].base[3];
			}
		};
		*out << std::endl;	
		site_=thissite;
	} else return UNEXPECTED;
	return 0;
};

/** @brief opens a .pro file in the modes "r" or "w".
  * @returns a pointer to the profile
**/
profile* profile::open(const char* filename, const char mode){
	sorted_[0]=0;
	sorted_[1]=1;
	sorted_[2]=2;
	sorted_[3]=3;
	sorted_[4]=4;
	open_=false;
	in=NULL;
	out=NULL;

	switch ( mode ) {
		case 'r':
			inFile.open(filename, std::fstream::in);
			//inFile.open(filename, std::fstream::in, std::ios::binary);
			if (!inFile.is_open() ){
				std::cerr << "cannot open " << filename << " for reading (1)." << std::endl;				
				exit(0);
			};
			in=&inFile;
			if (readheader()==BADHEADER){
				std::cerr << "cannot read header on " << filename << " (1). " << std::endl;				
				std::cerr << "Vesions of mapgd >=2.0 require headers on .pro files" << std::endl;				
			};
			break;
		case 'w':
			outFile.open(filename, std::ofstream::out);
			//outFile.open(filename, std::ofstream::out, std::ios::binary);
			if (!outFile.is_open() ){
				std::cerr << "cannot open " << filename << " for writing." << std::endl;				
				exit(0);
			};
			out=&outFile;
			break;
		default :
			std::cerr << "unkown filemode " << std::endl;
			exit(0);
	
	}
	open_=true;
	return this;
}

profile* profile::open(char mode)
{
	sorted_[0]=0;
	sorted_[1]=1;
	sorted_[2]=2;
	sorted_[3]=3;
	sorted_[4]=4;
	open_=false;
	
	switch ( mode ) {
		case 'r':
			in=&std::cin;
			if (readheader()==BADHEADER){
				std::cerr << "cannot find header in stdin (2). " << std::endl;				
				std::cerr << "Vesions of mapgd >=2.0 require headers on .pro files" << std::endl;				
			}
			break;
		case 'w':
			columns_=7;
			out=&std::cout;
			break;
		default :
			std::cerr << "unkown filemode " << std::endl;
			exit(0);
	
	};

	open_=true;
	donothing_=false;
	return this;
}
/** @brief closes a .pro file and unsets members.
  * @no return value
**/
void profile::close(void){
	open_=false;
};

profile::profile(){
	open_=false;
};

bool profile::is_open(void) const{
	return open_;
}

void profile::sort(void){
	count_t total[5]={0};
	for (unsigned int s=0; s<samples_;++s){
		if (site_.sample[s].masked) continue;
		total[0]+=site_.sample[s].base[0];
		total[1]+=site_.sample[s].base[1];
		total[2]+=site_.sample[s].base[2];
		total[3]+=site_.sample[s].base[3];
		total[4]+=site_.sample[s].base[4];
	};
	if (total[sorted_[0]]<total[sorted_[2]])
		std::swap(sorted_[0], sorted_[2]);
	if (total[sorted_[1]]<total[sorted_[3]])
		std::swap(sorted_[1], sorted_[3]);
	if (total[sorted_[2]]<total[sorted_[3]])
		std::swap(sorted_[2], sorted_[3]);
	if (total[sorted_[0]]<total[sorted_[1]])
		std::swap(sorted_[0], sorted_[1]);
	if (total[sorted_[1]]<total[sorted_[2]])
		std::swap(sorted_[1], sorted_[2]);
};

void profile::sort(count_t s){
	if (site_.sample[s].base[sorted_[0]]<site_.sample[s].base[sorted_[2]])
		std::swap(sorted_[0], sorted_[2]);
	if (site_.sample[s].base[sorted_[1]]<site_.sample[s].base[sorted_[3]])
		std::swap(sorted_[1], sorted_[3]);
	if (site_.sample[s].base[sorted_[2]]<site_.sample[s].base[sorted_[3]])
		std::swap(sorted_[2], sorted_[3]);
	if (site_.sample[s].base[sorted_[0]]<site_.sample[s].base[sorted_[1]])
		std::swap(sorted_[0], sorted_[1]);
	if (site_.sample[s].base[sorted_[1]]<site_.sample[s].base[sorted_[2]])
		std::swap(sorted_[1], sorted_[2]);
};
void profile::swap(count_t x, count_t y){	
	std::swap(sorted_[x], sorted_[y]);
};
count_t profile::getcount(count_t s, count_t c) const
{
	if (c<5) return site_.sample[s].base[sorted_[c]];
	else return 0;
};

const count_t *profile::getquartet(count_t s) const
{
	return site_.sample[s].base;
};

count_t profile::getcount(count_t c) const
{
	count_t total=0;
	if (c<5){
		for (unsigned int s=0; s<samples_;++s){
			if (site_.sample[s].masked) continue;
			total+=site_.sample[s].base[sorted_[c]];
		}
		return total;
	}
	else return 0;
};

count_t profile::getcoverage(count_t s) const
{
	if (s<samples_)	return site_.sample[s].base[0]+
				site_.sample[s].base[1]+
				site_.sample[s].base[2]+
				site_.sample[s].base[3];
	else {
		std::cerr << "Attempted to access a population that doesn't exist. Exiting." << std::endl;
		exit(0);
	};
};

count_t profile::getcoverage() const
{
	count_t total=0;
	for (unsigned int s=0; s<samples_;++s){
		if (site_.sample[s].masked) continue;
		total+=site_.sample[s].base[0]+
			site_.sample[s].base[1]+
			site_.sample[s].base[2]+
			site_.sample[s].base[3];
	};
	return total;
};

const count_t profile::getindex(count_t c) const
{
	if (c<5) return  sorted_[c];
	else return -1;
};

name_t profile::getname(count_t c) const
{
	if (c<5){
		return  profile::names_[sorted_[c]];
	}
	return '*';
};

name_t profile::getname_gt(count_t c) const
{
	if (c<5){
		if (getcount(c)==getcount(c+1) ) return '*';
		return  profile::names_[sorted_[c]];
	}
	return '*';
};

std::string profile::getids(void) const
{
	std::string str=site_.id0+'\t'+site_.id1;
	for (unsigned int x=0; x<site_.extra_ids.size(); ++x){
		str+='\t'+site_.extra_ids[x];
	};
	return str;
};

count_t profile::size(void) const
{
	return samples_;
};

count_t profile::maskedcount(void) const
{
	count_t count=0;
	for(std::vector<quartet>::const_iterator it = site_.sample.begin(); it != site_.sample.end(); ++it) if (it->masked) count++;
	return count;
};

void profile::maskall(void){
	for (unsigned int s=0; s<samples_;++s){
		site_.sample[s].masked=true;
	};
};

void profile::unmask(count_t a){
	site_.sample[a].masked=false;
};

void profile::unmask(quartet_t *q){
	q->masked=false;
}

void profile::mask(quartet_t *q){
	q->masked=true;
}

std::vector <quartet_t>::const_iterator profile::begin(void) const 
{
	return site_.sample.begin();
}	//returns a pointer to the end (unsorted)

std::vector <quartet_t>::iterator profile::begin(void)  
{
	return site_.sample.begin();
}	//returns a pointer to the end (unsorted)

std::vector <quartet_t>::const_iterator profile::end(void) const 
{
	return site_.sample.end();
}	//returns a pointer to the end (unsorted)

std::vector <quartet_t>::iterator profile::end(void)  
{
	return site_.sample.end();
}	//returns a pointer to the end (unsorted)

void profile::set_delim_quartet(const char &del)  
{
	delim_quartet=del;
}	

void profile::set_delim_column(const char &del)  
{
	delim_column=del;
}	


/*count_t major(const quartet_t a)
{
	?
	return 
};
count_t minor(const quartet_t b){
	?
	return 
};*/
count_t count(const quartet_t c){
	return c.base[0]+c.base[1]+c.base[2]+c.base[3];
}
