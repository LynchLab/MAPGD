#include "pro-file.h"

/*! \brief The default order of names_. Together with sorted_ member of Locus this specifies the identity of a quartet.
*/	
static const std::string names_="ACGTN";

/*! \brief The default mapping of quartet to names. memcpy this to reorder the sorted_ member of Locus. 
 */	
const count_t profile::defaultorder[5] = {0,1,2,3,4};

/*! \brief The default mapping of quartet to names.
 */	
const std::string profile::decodeid0(const id0_t &id){
	return header_.decodeid0(id);
}

/*! \brief takes a numerical id and returns the string in represents.
 */	
const std::string profile_header::decodeid0(const id0_t &id){
	if (lastid0==id) return lastid0_str;
	lastid0=id;
	if (id<id0.size() ) {
		lastid0_str=id0[id]; 
		return lastid0_str; 
	}
	else {
		std::cerr << __FILE__ << ":" << __LINE__ << ":" << id << ", " << lastid0 << ", " << id0.size() << ", " << id0[id] << "\n";
		exit(0);
		return "";
	}
//	return index.get_string(id);
}

/*! \brief takes a numerical id and returns the string in represents.
 */	
void profile::setsamples(count_t samples){
	header_.setsamples(samples);
}			

void profile::setcolumns(count_t a){
	header_.setcolumns(a);
}

const std::string profile::decodeid1(const id1_t &id){
	return header_.decodeid1(id);
}
const std::string profile_header::decodeid1(const id1_t &id){
	return std::to_string( (unsigned long long int)id );
}

const std::string profile::decodeextraid(const char &id, const size_t &a) {
	return header_.decodeextraid(id, a);
}
const std::string profile_header::decodeextraid(const char &id, const size_t &a) {
	return std::string(1, decodechar[size_t(id)]);
}

const id0_t profile::encodeid0(const std::string &id){
	return header_.encodeid0(id);
}

const id0_t profile_header::encodeid0(const std::string &id){
	if (lastid0_str==id) return lastid0;
	lastid0_str=id;
	std::map<std::string, count_t>::iterator search = id0_str.find(id);
	if(search != id0_str.end()) {
		lastid0=search->second;
		return lastid0;
	}
	else {
		control=(control|NEWID0);
		count_t emplace=id0_str.size();
		id0_str[id]=emplace;
		id0.push_back(id);
		lastid0=emplace;
		return lastid0;
	}
}

const id1_t profile::encodeid1(const std::string &id){
	return header_.encodeid1(id);
}

const id1_t profile_header::encodeid1(const std::string &id){
	return atoi(id.c_str() );
}

const char profile::encodeextraid(const char &id, const count_t &a){
	return header_.encodeextraid(id, a);
}
const char profile_header::encodeextraid(const char &id, const size_t &a){
	return encodechar[size_t(id)];
}

int profile::seek(std::streampos pos) {
	/*if (column==5 || in==NULL){
		std::cerr << "seek is currently only implemented on random access files opened in 6 or 7 column mode. Try re-running ei w/o the -s option.\n";
	};*/
	in->seekg(pos);
	read();
	return NONE;
}

void profile_header::clear(void){ column_names.clear(); }


profile_header::~profile_header(void){ column_names.clear(); }

int profile_header::setsamples(const count_t &samples) {
	*samples_=samples;
	column_names.clear();
	if (*columns_==6 or *columns_==7) column_names.push_back("scaffold");
	column_names.push_back("pos");
	if (*columns_==7) column_names.push_back("ref");

	for (unsigned int x=0; x<*samples_; ++x){
		column_names.push_back("sample_"+std::to_string( (unsigned long long int)(x+1) ) );
		sample_gof_.push_back(0);				// the number of samples (i.e. different individuals or populations) in the profile.
	};
	site_->resize(samples);
	return NONE;
}


int profile::copyheader(const profile &pro){
	header_=pro.header_;
	return NONE;
}

int profile_header::setcolumns(const count_t &x) {
	if (column_names.size()==0) setsamples(0);
	if (*columns_==5 and x==6){column_names.insert(column_names.begin(), "scaffold");}
	if (*columns_==5 and x==7){column_names.insert(column_names.begin(), "scaffold"); column_names.insert(column_names.begin()+2, "ref");}
	if (*columns_==6 and x==5){column_names.erase(column_names.begin() );}
	if (*columns_==6 and x==7){column_names.insert(column_names.begin()+2, "ref");}
	if (*columns_==7 and x==5){column_names.erase(column_names.begin() ); column_names.erase(column_names.begin()+1);}
	if (*columns_==7 and x==6){column_names.erase(column_names.begin()+2);}
	*columns_=x;
	return NONE;
}

const float_t profile::getsample_property(const count_t &a) const {
	return header_.getsample_property(a);
}
const float_t profile_header::getsample_property(const count_t &a) const {
	return sample_gof_[a];
}


const std::string profile_header::get_sample_name(const count_t &a) const{
	switch (*columns_){
		case 5:
			return get_column_name(a+1);
		break;
		case 6:
			return get_column_name(a+2);
		break;
		case 7:
			return get_column_name(a+3);
		break;
		default:
		break;
	};
	return NULL;
}

const std::string profile::get_sample_name(const count_t &a) const{
	return header_.get_sample_name(a);
};

int profile::set_sample_name(const std::string &str){
	switch (columns_){
		case 5:
			for (size_t a=0; a<samples_; ++a) if (!(header_.set_column_name(a+1, str) ) ) return BADHEADER;
			return 0;
		break;
		case 6:
			for (size_t a=0; a<samples_; ++a) if (!(header_.set_column_name(a+2, str) ) ) return BADHEADER;
			return 0;
		break;
		case 7:
			for (size_t a=0; a<samples_; ++a) if (!(header_.set_column_name(a+3, str) ) ) return BADHEADER;
			return 0;
		break;
	};
	return BADHEADER;
}

int profile::set_sample_name(const count_t &a, const std::string &str){
	switch (columns_){
		case 5:
			return header_.set_column_name(a+1, str);
		break;
		case 6:
			return header_.set_column_name(a+2, str);
		break;
		case 7:
			return header_.set_column_name(a+3, str);
		break;
	};
	return BADHEADER;
}

const std::string profile_header::get_column_name(const size_t &x) const{
	if (x<column_names.size() ) return column_names[x];
	std::cerr << __FILE__ << ":" << __LINE__ << ": attempted to access column that does not exist.\n";
	exit(0);
	return "ACCESS ERROR";
};

int profile_header::set_column_name(const size_t &x, const std::string &str) {
	if (x<column_names.size() ) {column_names[x]=str; return NONE;}
	return BADHEADER;
}

profile_header::profile_header(){}

void profile_header::init(profile *pro){
	sig_=&pro->sig_;			// were alleles thrown out if the allele only occurred in reads from one direction?
	read_=&pro->read_;			// file mode flag;
	write_=&pro->write_;			// file mode flag;
	binary_=&pro->binary_;			// file mode flag;
	mpileup_=&pro->mpileup_;
	noheader_=&pro->noheader_;
	delim_column=&pro->delim_column;	// the delimiter which seperates columns
	delim_quartet=&pro->delim_quartet;	// the delimiter that seperates counts in a quartet
	columns_=&pro->columns_;		// 5|6|7|more?
	samples_=&pro->samples_;		// the number of samples (i.e. different individuals or populations) in the profile.
	site_=&pro->site_;			// the number of samples (i.e. different individuals or populations) in the profile.
	size_=&pro->size_;			// the number of samples (i.e. different individuals or populations) in the profile.
	
	for (int x=0; x<256; x++) encodechar[x]=5;
	for (int x=0; x<256; x++) decodechar[x]='N';

	encodechar['A']=0;
	encodechar['a']=0;
	encodechar['C']=1;
	encodechar['c']=1;
	encodechar['G']=2;
	encodechar['g']=2;
	encodechar['T']=3;
	encodechar['t']=3;

	decodechar[0]='A';
	decodechar[1]='C';
	decodechar[2]='G';
	decodechar[3]='T';
	decodechar[4]='N';
}

profile_header::profile_header(profile *pro){
	init(pro);
}

int profile_header::writeheader(std::ostream *out){
	if (out==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempted writeheader when no ostream was open." << std::endl;
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
		*out << "@PR" << *delim_column << "VN:" << VER << *delim_column << "CN:" << *columns_ << *delim_column << "SN:" << *samples_ << *delim_column<< "MD:" << *binary_ <<std::endl;
		*out << "@ID";
		for (unsigned int x=0; x<column_names.size(); ++x)
			*out << *delim_column << column_names[x];
		*out << std::endl;
		notdone_=false;
	}
	return NONE;
}


int profile_header::writetailer(std::ostream *out){
	if (out==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": called writetailer when no ostream was open.\n";
		exit(UNEXPECTED);
	};
	bool notdone_=true;
	while(notdone_){
		if (*binary_){
			char e=EOBIN|control;
			out->write((char *)&(e), sizeof(char) );
			*out << '\n';
		}
		*out << "@LN:" << *size_ << std::endl;
		notdone_=false;
	}
	return NONE;
}

profile_header & profile_header::operator =(const profile_header& arg){
	if (this != &arg) { 
		*delim_column=*arg.delim_column;	// the delimiter which seperates columns
 		*delim_quartet=*arg.delim_quartet;	// the delimiter that seperates counts in a quartet
		*columns_=*arg.columns_;		// 5|6|7|more?
		*site_=*arg.site_;			// a vector to store the quartet.
		*samples_=*arg.samples_;		// the number of samples (i.e. different individuals or populations). This needs to be left to Locus.size().
		setsamples(*samples_);			// the right way to setsamples.
		column_names=arg.column_names;


		id0_str=arg.id0_str;
		id0=arg.id0;
//		index=arg.index;
		extraids=arg.extraids;			//extra ids associated with the quartet. (ref base identiy?).   
		sample_gof_=arg.sample_gof_;		// the number of samples (i.e. different individuals or populations) in the profile.

		lastid0=-1;				//initilize to 0-1;
	        lastid0_str="";				//initilize to "";

		memcpy (encodechar,arg.encodechar, 256*sizeof(count_t) );
		memcpy (decodechar,arg.decodechar, 256*sizeof(char) );
	}

	return *this;
};


//TODO FIX THIS, IT'S UGLY.

#define H_PARAM	1
#define H_ID	2
#define H_VER	3
#define H_COL	4
#define H_ERR	5
#define H_SAM	6
#define H_MOD	7
#define H_LN	8

int hash(std::string str){
	if (str=="@PR") return H_PARAM;
	if (str=="@ID") return H_ID;
	if (str=="VN") return H_VER;
	if (str=="CN") return H_COL;
	if (str=="SN") return H_SAM;
	if (str=="MD") return H_MOD;
	if (str=="@LN") return H_LN;
	else return H_ERR;
};

int profile::readheader(void){
	return header_.readheader(in);
}

int profile::writeheader(void)
{
	return header_.writeheader(out);
}

int profile_header::readheader(std::istream *in)
{
	//Headers should contain the following feilds:
	//The version number of proview used to create the file,
	//The number of columns in the file
	//The number of ids in the files
	//Header lines begin with the '@' character (a nod to samtools).
	if (in==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": called readheader when no istream was open.\n";
		exit(UNEXPECTED);
	};
	std::string line;
	std::vector <std::string> column, args, arg;

	control=0;

	bool notdone_=true, reading_id=false, reading_pr=true;

#ifdef DEBUG 
	std::cerr << "attempting to read header.\n";
#endif

	while(notdone_){
		//"VN" version number
		//"@IDs" Names of fields.
		if (std::getline(*in, line) ){
			if ( line[0]!='@' || *noheader_ ){
#ifdef DEBUG 
				std::cerr << "bad header format.\n";
#endif
				column=split(line, *delim_column);
				if (column[1][0]=='s') {
#ifdef DEBUG 
					std::cerr << "assuming Takahiro type file.\n"; 
#endif
					*noheader_=1;
}
				if ( *noheader_ ) {
					*delim_quartet='/';
					*columns_=7; 
					switch (*columns_){
						case 5:
							setsamples(column.size()-1);
						break;
						case 6:
							setsamples(column.size()-2);
						break;
						case 7:
							site_->set_extraid(5, 0);
							setsamples(column.size()-3);
						break;
					};
					for (unsigned int x=0; x<column.size(); ++x) column_names[x]=column[x];
					notdone_=false;
				}
				else {
					*mpileup_=true;
					*samples_=(column.size()-3)/3;
					setsamples(*samples_);
					in->putback('\n');
					site_->set_extraid(5, 0);
					*columns_=7;
					for (std::string::reverse_iterator rit=line.rbegin(); rit!=line.rend(); ++rit) in->putback(*rit);
					notdone_=false;
				}
			} else {
				args=split(line, *delim_column);
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
							*columns_=atoi(arg[1].c_str() );
							if (*columns_==7) site_->set_extraid(5, 0);
							break;
						case H_SAM:
							*samples_=atoi(arg[1].c_str() );
							setsamples(*samples_);
							break;
						case H_MOD:
							*binary_=atoi(arg[1].c_str() );
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
								std::cerr << __FILE__ << ":" << __LINE__ << ": unexpected field encountered in header (" << *argit << "). File may not have opened correctly. Try removing this field and running the program again." << std::endl;
								return BADHEADER;
							}
						break;
					}
				}
				if (in!=&std::cin){
					std::streampos S=in->tellg();
					in->seekg(0, in->end);
					in->clear();
					in->unget();
					in->unget();
					while( (in->peek() )!='\n') {in->unget();}
					if ( not (std::getline(*in, line) ) ) { 
						std::cerr << __FILE__ << ":" << __LINE__ << ": encountred unexpected EOF. Is the file truncated? (1)" << std::endl;
						return BADHEADER; 
					}
					if ( not (std::getline(*in, line) ) ) { 
						std::cerr << __FILE__ << ":" << __LINE__ << ": encountred unexpected EOF. Is the file truncated? (2)" << std::endl;
						return BADHEADER; 
					}
					args=split(line, *delim_column);
					if (args.size()==0) { 
						std::cerr << __FILE__ << ":" << __LINE__ << ": line miss-formated. Is the file truncated?" << std::endl; 
						return BADHEADER;
					}
					for(std::vector<std::string>::iterator argit = args.begin(); argit != args.end(); ++argit){
						arg=split(*argit, ':');
						switch(hash(arg[0]) ){
							case H_LN:
								*size_=atoi(arg[1].c_str() );
							break;
							default:
							break;
						}
					}
					in->seekg(S);
					in->clear();
				};
			} 
		} else {
			std::cerr << __FILE__ << ":" << __LINE__ << ": encountred unexpected EOF. Is the file empty?" << std::endl;
			return BADHEADER;
		}; 
	}
#ifdef DEBUG 
	std::cerr << "done reading header.\n";
#endif
	return NONE;
}

int profile::read()
{
	return read(site_);
}

int profile::read(Locus &site)
{
	while (site.sample.size()<samples_) site.sample.push_back(quartet_t() );
	while (site.extraid.size()<site_.extraid.size() ) site.extraid.push_back(0);
	if(binary_) { 
		return readb(site);
	} else {
		if (mpileup_) return readm(site);
		else return readt(site);
	};
}

int profile::readm(Locus &site) 
{
	std::string line;
	memcpy(site.sorted_, defaultorder, 5*sizeof(count_t) );
	std::vector <std::string> column, quartet;
	

	if( std::getline(*in, line) ){
		column = split(line, delim_column);
		if(column.size() < 5){
			std::cerr << "WARNING:" << __FILE__ << ":" << __LINE__ << ": Skipping line " << getlinenumber() << " with only " << column.size() << " fields." << std::endl;
			return 0;
		}

		std::vector <quartet_t>::iterator it=site.sample.begin();
		std::vector <quartet_t>::iterator it_end=site.sample.end();
		std::vector <std::string>::iterator column_it=column.begin();	

		site.id0=encodeid0(column[0]);
		site.id1=encodeid1(column[1]);
		site.extraid[0]=encodeextraid(column[2][0], 0);

		column_it+=4;
		while (it!=it_end){
    			memset(it->base, 0, sizeof(count_t)*5 );
			scan(site, *column_it, *it);
			column_it+=3;
			++it;
		}
		return 0;
		
	} else {
		return EOF;
	}
	return UNEXPECTED;
}

void inline profile::scan(const Locus & site, const std::string &str, quartet_t &q)
{
	size_t skip;
	std::string number;
	std::string::const_iterator it=str.begin();
	std::string::const_iterator end=str.end();
	size_t lim=0;
	while (it!=end) {
		if (*it == '^') {it++; if(it==end) break; it++; continue;}
		if (*it == '$' || *it=='*') {it++; continue;}
		else if (*it == '+' || *it == '-') {
			number = "";
			it++;
			while(isdigit(*it) ) {
				number.push_back(*it); 
				it++; 
			}
			lim=atoi(number.c_str() );
 			for(skip=0; skip<lim; ++skip) it++;
			continue;
    		} else if(*it == ',' || *it=='.') {
			q.base[site.extraid[0] ]++;
			#ifdef DEBUG
				if (site.extraid[0]==5) {
					std::cerr << __FILE__ << ":" << __LINE__ << ": no reference nucleotide set.\n";
				}
			#endif 
		} else {
			q.base[size_t(header_.encodeextraid( (char) (*it), 0) )]++;
		}
		it++;
	}
}

int profile::readb(Locus &site){
	memcpy(site.sorted_, defaultorder, 5*sizeof(count_t) );
	in->read((char *)&(header_.control), sizeof(char) );
	if (header_.control){
		if (header_.control&NEWID0){
			std::string str;
			std::getline(*in, str);
			encodeid0(str);
		}
		if (header_.control&EOBIN) return EOF;
	}
	switch (columns_){
		case 5:
		case 6:
			in->read ( (char *)&site.id0, sizeof(id0_t) );
			in->read ( (char *)&site.id1, sizeof(id1_t) );
			for (unsigned int x=0; x<samples_; ++x) in->read ( (char *)site.sample[x].base, 4*sizeof(count_t) );
		break;
		case 7:
			in->read ( (char *)&site.id0, sizeof(id0_t) );
			in->read ( (char *)&site.id1, sizeof(id1_t) );
			in->read ( (char *)&site.extraid[0], sizeof(char) );
			for (unsigned int x=0; x<samples_; ++x) in->read ( (char *)site.sample[x].base, 4*sizeof(count_t) );
		break;
		default:
			std::cerr << __FILE__ << ":" << __LINE__ << ": error in binary file formating. "<< std::endl;
			exit(0);
		break;	
	}
	if (in->eof() ) return EOF;
	else return 0;
};

int profile::readt(Locus &site){
	//The order of nucleotides read from a quartet file. Replace this with a memcopy of a constant.
	std::string line;
	memcpy(site.sorted_, defaultorder, 5*sizeof(count_t) );
	std::vector <std::string> column, quartet;

	if (in==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempted read line when no istream was open." << std::endl;
		exit(UNEXPECTED);
	};

	if (std::getline(*in, line) ){
		if (line[0]=='@') return EOF;
		switch (columns_){
		//READ A PRO FILE
			case 5:
			/* check for start of new scaffold*/
				if (line[0]!='>'){
					column=split(line, delim_column);
					if ( (column.size()-1) != samples_){
						std::cerr << __FILE__ << ":" << __LINE__ << ": could not parse line : \"" << line << "\"" << std::endl;
						std::cerr << column.size() << " fields found." << std::endl;
						std::cerr << "delimiter : \'" << delim_column << "\'" << std::endl;
						std::cerr << "columns : \'" << columns_ << "\'" << std::endl;
						exit(0);
						return EOF;
					}
					site.id1=encodeid1(column[0]);
					for (unsigned int x=0; x< samples_; ++x){
						quartet=split(column[x+1], delim_quartet );
						site.sample[x].base[0]=atoi(quartet[0].c_str() );
						site.sample[x].base[1]=atoi(quartet[1].c_str() );
						site.sample[x].base[2]=atoi(quartet[2].c_str() );
						site.sample[x].base[3]=atoi(quartet[3].c_str() );
					};
				} else {
					line.erase(line.begin());
					site.id0=encodeid0(line);
					read(site);
				}
				
			break;
			case 6:
			column=split(line, delim_column);
			if ( (column.size()-2) != samples_){
				std::cerr << __FILE__ << ":" << __LINE__ << ": could not parse line : \"" << line << "\"" << std::endl;
				std::cerr << column.size() << " fields found." << std::endl;
				std::cerr << "delimiter : \'" << delim_column << "\'" << std::endl;
				std::cerr << "columns : \'" << columns_ << "\'" << std::endl;
				exit(0);
				return EOF;
			};
			site.id0=encodeid0(column[0]);
			site.id1=encodeid1(column[1]);
			for (size_t x=0; x < samples_; ++x){
				quartet=split(column[x+2], delim_quartet);
				site.sample[x].base[0]=atoi(quartet[0].c_str() );
				site.sample[x].base[1]=atoi(quartet[1].c_str() );
				site.sample[x].base[2]=atoi(quartet[2].c_str() );
				site.sample[x].base[3]=atoi(quartet[3].c_str() );
			};
			break;
			case 7:
			column=split(line, delim_column);
			if ( (column.size()-3)!=samples_){
				std::cerr << __FILE__ << ":" << __LINE__ << ": could not parse line : \"" << line << "\"" << std::endl;
				std::cerr << column.size() << " fields found." << std::endl;
				std::cerr << "delimiter : \'" << delim_column << "\'" << std::endl;
				std::cerr << "columns : \'" << columns_ << "\'" << std::endl;
				exit(0);
				return EOF;
			};
			site.id0=encodeid0(column[0]);
			site.id1=encodeid1(column[1]);
			site.extraid[0]=encodeextraid(column[2][0], 0);
			for (size_t x=0; x<samples_; ++x){
				quartet=split(column[x+3], delim_quartet);
				site.sample[x].base[0]=atoi(quartet[0].c_str() );
				site.sample[x].base[1]=atoi(quartet[1].c_str() );
				site.sample[x].base[2]=atoi(quartet[2].c_str() );
				site.sample[x].base[3]=atoi(quartet[3].c_str() );
			}
			break;
		}
		return 0;
	}
	return EOF;
}

int profile::copy(const profile &pro){
	site_=pro.site_;
	return 0;
}

int profile::write(void){
	return write(site_);	
}

int profile::write(const Locus &site){
	size_+=1;
	if (binary_) return writeb(site);
	else return writet(site);
}

int profile::writeb (const Locus &thissite){
	if (out==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempted write when no ostream was open." << std::endl;
		exit(UNEXPECTED);
	}
	out->write((char *)&(header_.control), sizeof(char) );
	if ( (header_.control)&NEWID0) *out << decodeid0(thissite.id0) << std::endl;
        switch (columns_){
                case 5:
                case 6:
                        out->write((char *)&thissite.id0, sizeof(id0_t) );
			out->write((char *)&thissite.id1, sizeof(id1_t) );
			for (unsigned int x=0; x<samples_; ++x) out->write((char *)thissite.sample[x].base, 4*sizeof(count_t) );
                break;
                case 7:
                        out->write((char *)&thissite.id0, sizeof(id0_t) );
			out->write((char *)&thissite.id1, sizeof(id1_t) );
			out->write((char *)thissite.extraid[0], sizeof(char) );
			for (unsigned int x=0; x<samples_; ++x) out->write((char *)thissite.sample[x].base, 4*sizeof(count_t) );
                break;
                default:
			std::cerr << __FILE__ << ":" << __LINE__ << ": error in binary file formating. "<< std::endl;
                        exit(0);
                break;
        }
	(header_.control)=0;
	site_=thissite;
        return 0;
}

int profile::writet(const Locus &thissite){
	if (out==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempted write when no outstream was open." << std::endl;
		exit(UNEXPECTED);
	};
	switch (columns_){
	case 5:
		if (thissite.id0!=site_.id0){
			*out << '>' << decodeid0(thissite.id0) << std::endl;
		};
		*out << decodeid1(thissite.id1);
		for (size_t x=0; x < thissite.sample.size(); ++x){
			*out << delim_column;
			*out << thissite.sample[x];
		};
		*out << std::endl;
		site_=thissite;
	break;
	case 6:
		*out << decodeid0(thissite.id0) << delim_column << decodeid1(thissite.id1);
		for (size_t x=0; x < thissite.sample.size(); ++x){
			*out << delim_column;
			*out << thissite.sample[x];
		};
		*out << std::endl;	
		site_=thissite;
	break;
	case 7: 
		*out << decodeid0(thissite.id0) << delim_column << decodeid1(thissite.id1) << delim_column << decodeextraid(thissite.get_extraid(0), 0);
		for (size_t x=0; x < thissite.sample.size(); ++x){
			*out << delim_column;
			*out << thissite.sample[x];
		};
		*out << std::endl;	
		site_=thissite;
	break;
	default:
		return UNEXPECTED;
	break;
	}
	return 0;
};

/** @brief opens a .pro file in the modes "r" or "w".
  *
**/
void profile::open(const char* filename, std::ios_base::openmode mode){

	memcpy(site_.sorted_, defaultorder, 5*sizeof(count_t) );
	open_=false;
	in=NULL;
	out=NULL;

	if ( mode & std::fstream::in ){
		inFile.open(filename, std::ifstream::in);
		if (!inFile.is_open() ){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open " << filename << " for reading (1)." << std::endl;				
			exit(0);
		};
		in=&inFile;
		if (readheader()==BADHEADER){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot read header on " << filename << " (1). " << std::endl;
			std::cerr << "Vesions of mapgd >=0.3 require headers on .pro files" << std::endl;
			exit(0);
		};
		read_=true;
	} else if ( mode & std::fstream::out ){
		outFile.open(filename, std::ofstream::out);
		if (!outFile.is_open() ){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open " << filename << " for writing." << std::endl;				
			exit(0);
		};
		out=&outFile;
		write_=true;
	};
	if (mode & std::fstream::binary) binary_=true;
	open_=true;
}

void profile::open(std::ios_base::openmode mode)
{
	open_=false;
	if ( mode & std::fstream::in ){
		in=&std::cin;
		if (readheader()==BADHEADER){
			std::cerr << "cannot find header in stdin (2). " << std::endl;				
			std::cerr << "Vesions of mapgd >=0.3 require headers on .pro files" << std::endl;			
			exit(0);
		}
		read_=true;
	} else if ( mode & std::fstream::out) {
		out=&std::cout;
		write_=true;
	}
	if (mode & std::fstream::binary) binary_=true;
	open_=true;
	donothing_=false;
}

/** @brief closes a .pro file and unsets members.
  * @no return value
**/

void profile::close(void){
	if (write_){
		if (not noheader_) header_.writetailer(out);
		if (outFile.is_open() ) outFile.close();
	}
	if(read_) if (inFile.is_open() ) inFile.close();
	header_.clear();
	read_=false;
	write_=false;
	open_=false;
};

profile::profile(){
	open_=false;
	read_=false;
	write_=false;
	binary_=false;
	mpileup_=false;
	noheader_=false;

	out=NULL;
	in=NULL;

	size_=0;
	samples_=0;
	columns_=7;
	delim_column='\t';
	delim_quartet='/';
	
	header_.init(this);
};

bool profile::is_open(void) const{
	return open_;
}

void profile::sort(void){
	site_.sort();
}

std::string profile::getids(void)
{
	std::string str=decodeid0(site_.id0)+'\t'+decodeid1(site_.id1);
	for (unsigned int x=0; x<site_.extraid.size(); ++x){
		str+='\t'+decodeextraid(site_.extraid[x], x);
	};
	return str;
};

std::string profile::getids(const Locus &site)
{
	std::string str=decodeid0(site.id0)+'\t'+decodeid1(site.id1);
	for (unsigned int x=0; x<site.extraid.size(); ++x){
		str+='\t'+decodeextraid(site.extraid[x], x);
	};
	return str;
};
	
size_t profile::size(void) const
{
	return samples_;
};

void profile::maskall(void){
	site_.maskall();
};

void profile::unmaskall(void){
	site_.maskall();
};

count_t profile::maskedcount(void){
	return site_.maskedcount();
};

/*
*/

void profile::set_delim_quartet(const char &del)  
{
	delim_quartet=del;
}	

void profile::set_delim_column(const char &del)  
{
	delim_column=del;
}	

void profile::setbasecount(const count_t &x, const count_t &y, const count_t &z)
{
	site_.sample[x].base[y]=z;
}

const id1_t profile::getlinenumber(void) const
{
	return size_;
}

const id0_t profile::get_id0(void) const
{
	return site_.id0;
}

const id1_t profile::get_id1(void) const 
{
	return site_.id1;
}

void profile::set_id0(const id0_t &id0) 
{
	site_.id0=id0;
}

void profile::set_id1(const id1_t &id1) 
{
	site_.id1=id1;
}

const char profile::get_extraid(const size_t &x) const 
{
	return site_.get_extraid(x);
}

void profile::set_extraid(const char &eid, const size_t &x)
{
	site_.set_extraid(eid, x);
}

//===WHOOPS SOMETHING IS WRONG HERE!
	
char Locus::getname(const count_t &c) const
{
        if (c<5){
                if (getcount(c)==0 ) return '.';
                return  names_[sorted_[c]];
        }
        return '.';
}

char Locus::getname_gt(const count_t &c) const
{
        if (c<5){
                if (getcount(c)==getcount(c+1) || getcount(c)==0 ) return '.';
                return  names_[sorted_[c]];
        }
        return '.';
}
