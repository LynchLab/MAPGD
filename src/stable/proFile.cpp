#include "proFile.h"

/*Design notes: This library is probably poorly thought out. More or less everything works through the 'profile'. However, 
  It may make a whole lot more sense for the 'profile' to be a file like structure that just co-ordinates reading, writing, 
  formating and all that jazz. Most of the work should probably be done by the "site" structure.
 */

/** @brief Read a .pro file in the five column format.
  * @returns 0 iff successful
**/

const std::string profile::names_ ="ACGT*";

const count_t profile::defaultorder[5] = {0,1,2,3,4};

site_t::site_t(void){
	sample.clear();
	samples_=0;
	id0=0;
	id1=0;	
};

site_t::~site_t(void){
	if (sample.size()>0) sample.clear();
};


site_t::site_t(count_t size){
	sample.assign(size, quartet() );
}

site_t & site_t::operator =(const site_t& arg){
        sample=arg.sample; 	
        id0=arg.id0;
	id1=arg.id1;  
        extraid=arg.extraid;   
}

const std::string profile::decodeid0(const count_t &id){
	return header_.decodeid0(id);
}
const std::string profile_header::decodeid0(const count_t &id){
	if (lastid0==id) return lastid0_str;
	id==lastid0; 
	lastid0_str=id0[id];
	return lastid0_str;
}

int profile::setsamples(count_t samples){
	return header_.setsamples(samples);
}			


int profile::setcolumns(count_t a){
	return header_.setcolumns(a);
}
const std::string profile::decodeid1(const uint64_t &id){
	return header_.decodeid1(id);
}
const std::string profile_header::decodeid1(const uint64_t &id){
	return std::to_string(id);
}

const std::string profile::decodeextraid(const count_t &id, const count_t &a) {
	return header_.decodeextraid(id, a);
}
const std::string profile_header::decodeextraid(const count_t &id, const count_t &a) {
	return decodechar[id];
}

const count_t profile::encodeid0(const std::string &id){
	return header_.encodeid0(id);
}

const count_t profile_header::encodeid0(const std::string &id){
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

const uint64_t profile::encodeid1(const std::string &id){
	return header_.encodeid1(id);
}

const uint64_t profile_header::encodeid1(const std::string &id){
	return atoi(id.c_str() );
}

const count_t profile::encodeextraid(const char &id, const count_t &a){
	return header_.encodeextraid(id, a);
}
const count_t profile_header::encodeextraid(const char &id, const count_t &a){
	return encodechar[id];
}

int profile::seek(std::streampos pos) {
	/*if (column==5 || in==NULL){
		std::cerr << "seek is currently only implemented on random access files opened in 6 or 7 column mode. Try re-running ei w/o the -s option.\n";
	};*/
	in->seekg(pos);
	read();
	return NONE;
}

profile_header::~profile_header(void){
	column_names.clear();
};

int profile_header::setsamples(const count_t &samples) {
	*samples_=samples;
	column_names.clear();
	if (*columns_==6 or *columns_==7) column_names.push_back("scaffold");
	column_names.push_back("pos");
	if (*columns_==7) column_names.push_back("ref");

	for (unsigned int x=0; x<*samples_; ++x){
		column_names.push_back("sample_"+std::to_string( (unsigned long long int)(x+1) ) );
		site_->sample.push_back(quartet_t() );
		sample_gof_.push_back(0);				// the number of samples (i.e. different individuals or populations) in the profile.
	};
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


const std::string profile_header::getsample_name(const count_t &a) const{
	switch (*columns_){
		case 5:
			return getcolumn_name(a+1);
		break;
		case 6:
			return getcolumn_name(a+2);
		break;
		case 7:
			return getcolumn_name(a+3);
		break;
	};
}

const std::string profile::getsample_name(const count_t &a) const{
	return header_.getsample_name(a);
};

int profile::setsample_name(const count_t &a, const std::string &str){
	switch (columns_){
		case 5:
			header_.setcolumn_name(a+1, str);
		break;
		case 6:
			header_.setcolumn_name(a+2, str);
		break;
		case 7:
			header_.setcolumn_name(a+3, str);
		break;
	};
}

const std::string profile_header::getcolumn_name(const count_t &x) const{
	return column_names[x];
};

int profile_header::setcolumn_name(const count_t &x, const std::string &str) {
	column_names[x]=str;
	return NONE;
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
		std::cerr << "attempted writetailer when no outstream was open." << std::endl;
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
        //*sig_=*arg.sig;                             // were alleles thrown out if the allele only occurred in reads from one direction?

        *delim_column=*arg.delim_column;                             // the delimiter which seperates columns
        *delim_quartet=*arg.delim_quartet;                            // the delimiter that seperates counts in a quartet
        *columns_=*arg.columns_;                         // 5|6|7|more?
        site_t *site_;                                  // a vector to store the calls from reads
        *samples_=*arg.samples_;                         // the number of samples (i.e. different individuals or populations) in the profile.

	*columns_=*arg.columns_;
	setsamples(*samples_);
	column_names.clear();

	for (unsigned int x=0; x<arg.column_names.size(); ++x){
		column_names.push_back(arg.column_names[x]);
	};
	//TODO Need to copy over maps and what not too!!	

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
	header_.readheader(in);
}

int profile::writeheader(void)
{
	header_.writeheader(out);
}

int profile_header::readheader(std::istream *in)
{
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

	control=0;

	bool notdone_=true, reading_pr=false, reading_id=false;

#ifdef DEBUG 
	std::cerr << "attempting to read header.\n";
#endif

	while(notdone_){
		//"VN" version number
		//"@IDs" Names of fields.
		if (std::getline(*in, line)!=NULL){
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
							site_->extraid.push_back(0);
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
					site_->extraid.push_back(0);
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
							if (*columns_==7) site_->extraid.push_back(0);
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
								std::cerr << "Warning: unexpected field encountered in header (" << *argit << "). File may not have opened correctly. Try removing this field and running the program again." << std::endl;
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
					std::getline(*in, line);
					std::getline(*in, line);
					args=split(line, *delim_column);
					if (args.size()==0) return BADHEADER;
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
			std::cerr << "Error: encountred unexpected EOF. Is the file empty?" << std::endl;
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

int profile::read(site_t &site)
{
	while (site.sample.size()<samples_) site.sample.push_back(quartet_t() );
	while (site.extraid.size()<site_.extraid.size() ) site.extraid.push_back(0);
	site.samples_=samples_;
	if(binary_) { 
		return readb(site);
	} else {
		if (mpileup_) return readm(site);
		else return readt(site);
	};
}

int profile::readm(site_t &site) 
{
	std::string line;
	memcpy(site.sorted_, defaultorder, 5*sizeof(count_t) );
	std::vector <std::string> column, quartet;
	

	if( std::getline(*in, line)!=NULL){
		column = split(line, delim_column);
		if(column.size() < 5){
			std::cerr << "WARNING[proFile.cpp]: Skipping line " << getlinenumber() << " with only " << column.size() << " fields." << std::endl;
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

void inline profile::scan(const site_t & site, const std::string &str, quartet_t &q)
{
	count_t j;
	std::string number;
	std::string::const_iterator it=str.begin();
	std::string::const_iterator end=str.end();
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
 			for(j=0; j<atoi(number.c_str() )-1; ++j) it++;
			continue;
    		} else if(*it == ',' || *it=='.') {
			q.base[site.extraid[0] ]++;
		} else {
			q.base[header_.encodeextraid( (count_t)*it, 0) ]++;
		}
		it++;
	}
}

int profile::readb(site_t &site){
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
			in->read ( (char *)&site.id0, sizeof(count_t) );
			in->read ( (char *)&site.id1, sizeof(uint64_t) );
			for (unsigned int x=0; x<samples_; ++x) in->read ( (char *)site.sample[x].base, 4*sizeof(count_t) );
		break;
		case 7:
			in->read ( (char *)&site.id0, sizeof(count_t) );
			in->read ( (char *)&site.id1, sizeof(uint64_t) );
			in->read ( (char *)&site.extraid[0], sizeof(count_t) );
			for (unsigned int x=0; x<samples_; ++x) in->read ( (char *)site.sample[x].base, 4*sizeof(count_t) );
		break;
		default:
			std::cerr << "Error in binary file formating. "<< std::endl;
			exit(0);
		break;	
	}
	if (in->eof() ) return EOF;
	else return 0;
};

int profile::readt(site_t &site){
	//The order of nucleotides read from a quartet file. Replace this with a memcopy of a constant.
	std::string line;
	memcpy(site.sorted_, defaultorder, 5*sizeof(count_t) );
	std::vector <std::string> column, quartet;

	if (in==NULL){
		std::cerr << "attempted read when no instream was open." << std::endl;
		exit(UNEXPECTED);
	};

	if (std::getline(*in, line)!=NULL){
		if (line[0]=='@') return EOF;
		switch (columns_){
		//READ A PRO FILE
			case 5:
			/* check for start of new scaffold*/
				if (line[0]!='>'){
					column=split(line, delim_column);
					if ( (column.size()-1) != samples_){
						std::cerr << "could not parse line : \"" << line << "\"" << std::endl;
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
				std::cerr << "could not parse line : \"" << line << "\"" << std::endl;
				std::cerr << column.size() << " fields found." << std::endl;
				std::cerr << "delimiter : \'" << delim_column << "\'" << std::endl;
				std::cerr << "columns : \'" << columns_ << "\'" << std::endl;
				exit(0);
				return EOF;
			};
			site.id0=encodeid0(column[0]);
			site.id1=encodeid1(column[1]);
			for (unsigned int x=0; x < samples_; ++x){
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
				std::cerr << "could not parse line : \"" << line << "\"" << std::endl;
				std::cerr << column.size() << " fields found." << std::endl;
				std::cerr << "delimiter : \'" << delim_column << "\'" << std::endl;
				std::cerr << "columns : \'" << columns_ << "\'" << std::endl;
				exit(0);
				return EOF;
			};
			site.id0=encodeid0(column[0]);
			site.id1=encodeid1(column[1]);
			site.extraid[0]=encodeextraid(column[2][0], 0);
			for (unsigned int x=0; x<samples_; ++x){
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
	write(site_);	
}

int profile::write(const site_t &site){
	size_+=1;
	if (binary_) return writeb(site);
	else return writet(site);
}
int profile::writeb (const site_t &thissite){
	if (out==NULL){
		std::cerr << "attempted write when no outstream was open." << std::endl;
		exit(UNEXPECTED);
	}
	out->write((char *)&(header_.control), sizeof(char) );
	if ( (header_.control)&NEWID0) *out << decodeid0(thissite.id0) << std::endl;
        switch (columns_){
                case 5:
                case 6:
                        out->write((char *)&thissite.id0, sizeof(count_t) );
			out->write((char *)&thissite.id1, sizeof(uint64_t) );
			for (unsigned int x=0; x<samples_; ++x) out->write((char *)thissite.sample[x].base, 4*sizeof(count_t) );
                break;
                case 7:
                        out->write((char *)&thissite.id0, sizeof(count_t) );
			out->write((char *)&thissite.id1, sizeof(uint64_t) );
			out->write((char *)thissite.extraid[0], sizeof(count_t) );
			for (unsigned int x=0; x<samples_; ++x) out->write((char *)thissite.sample[x].base, 4*sizeof(count_t) );
                break;
                default:
                        std::cerr << "Error in binary file formating. "<< std::endl;
                        exit(0);
                break;
        }
	(header_.control)=0;
	site_=thissite;
        return 0;
}

int profile::writet(const site_t &thissite){
	if (out==NULL){
		std::cerr << "attempted write when no outstream was open." << std::endl;
		exit(UNEXPECTED);
	};
	switch (columns_){
	case 5:
		if (thissite.id0!=site_.id0){
			*out << '>' << decodeid0(thissite.id0) << std::endl;
		};
		*out << decodeid1(thissite.id1);
		for (unsigned int x=0; x < samples_; ++x){
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
	break;
	case 6:
		*out << decodeid0(thissite.id0) << delim_column << decodeid1(thissite.id1);
		for (unsigned int x=0; x < samples_; ++x){
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
	break;
	case 7: 
		*out << decodeid0(thissite.id0) << delim_column << decodeid1(thissite.id1) << delim_column << decodeextraid(thissite.extraid[0], 0);
		for (unsigned int x=0; x < samples_; ++x){
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
	break;
	default:
		return UNEXPECTED;
	break;
	}
	return 0;
};

/** @brief opens a .pro file in the modes "r" or "w".
  * @returns a pointer to the profile
**/
profile* profile::open(const char* filename, const char *mode){

	memcpy(site_.sorted_, defaultorder, 5*sizeof(count_t) );
	open_=false;
	in=NULL;
	out=NULL;

	if ( mode=="r" ){
		inFile.open(filename, std::fstream::in);
		if (!inFile.is_open() ){
			std::cerr << "cannot open " << filename << " for reading (1)." << std::endl;				
			exit(0);
		};
		in=&inFile;
		if (readheader()==BADHEADER){
			std::cerr << "cannot read header on " << filename << " (1). " << std::endl;				
			std::cerr << "Vesions of mapgd >=2.0 require headers on .pro files" << std::endl;				
		};
		read_=true;
	} else if (mode=="w"){
		outFile.open(filename, std::ofstream::out);
		if (!outFile.is_open() ){
			std::cerr << "cannot open " << filename << " for writing." << std::endl;				
			exit(0);
		};
		out=&outFile;
		write_=true;
	} else if (mode=="wb"){
		outFile.open(filename, std::ofstream::out);
		if (!outFile.is_open() ){
			std::cerr << "cannot open " << filename << " for writing." << std::endl;				
			exit(0);
		};
		out=&outFile;
		write_=true;
		binary_=true;
	} else	{
		std::cerr << "unkown filemode " << std::endl;
		exit(0);
	}
	open_=true;
	return this;
}

profile* profile::open(const char *mode)
{
	open_=false;

	if (mode=="r"){
		in=&std::cin;
		if (readheader()==BADHEADER){
			std::cerr << "cannot find header in stdin (2). " << std::endl;				
			std::cerr << "Vesions of mapgd >=2.0 require headers on .pro files" << std::endl;				
		}
		read_=true;
	} else if (mode=="w") {
		out=&std::cout;
		write_=true;
	} else if (mode=="wb") {
		out=&std::cout;
		binary_=true;
		write_=true;
	} else{
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
	if (write_){
		header_.writetailer(out);
		if (outFile.is_open() ) outFile.close();
	}
	if(read_) if (inFile.is_open() ) inFile.close();
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
	columns_=5;
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
void site_t::sort(void){
	count_t total[5]={0};
	for (unsigned int s=0; s<sample.size();++s){
		if (sample[s].masked) continue;
		total[0]+=sample[s].base[0];
		total[1]+=sample[s].base[1];
		total[2]+=sample[s].base[2];
		total[3]+=sample[s].base[3];
		total[4]+=sample[s].base[4];
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

void site_t::sort(count_t s){
	if (sample[s].base[sorted_[0]]<sample[s].base[sorted_[2]])
		std::swap(sorted_[0], sorted_[2]);
	if (sample[s].base[sorted_[1]]<sample[s].base[sorted_[3]])
		std::swap(sorted_[1], sorted_[3]);
	if (sample[s].base[sorted_[2]]<sample[s].base[sorted_[3]])
		std::swap(sorted_[2], sorted_[3]);
	if (sample[s].base[sorted_[0]]<sample[s].base[sorted_[1]])
		std::swap(sorted_[0], sorted_[1]);
	if (sample[s].base[sorted_[1]]<sample[s].base[sorted_[2]])
		std::swap(sorted_[1], sorted_[2]);
};
void site_t::swap(count_t x, count_t y){	
	std::swap(sorted_[x], sorted_[y]);
};
count_t site_t::getcount(count_t s, count_t c) const
{
	if (c<5) return sample[s].base[sorted_[c]];
	else return 0;
};

const count_t * site_t::getquartet(count_t s) const
{
	return sample[s].base;
};

count_t site_t::getcount(count_t c) const
{
	count_t total=0;
	if (c<5){
		for (unsigned int s=0; s<samples_;++s){
			if (sample[s].masked) continue;
			total+=sample[s].base[sorted_[c]];
		}
		return total;
	}
	else return 0;
};

count_t site_t::getcoverage(count_t s) const
{
	if (s<samples_) return sample[s].base[0]+
				sample[s].base[1]+
				sample[s].base[2]+
				sample[s].base[3];
	else {
		std::cerr << "Attempted to access a population that doesn't exist. Exiting." << std::endl;
		exit(0);
	};
};

count_t site_t::getcoverage() const
{
	count_t total=0;
	for (unsigned int s=0; s<samples_;++s){
		if (sample[s].masked) continue;
		total+=sample[s].base[0]+
			sample[s].base[1]+
			sample[s].base[2]+
			sample[s].base[3];
	};
	return total;
};

const count_t site_t::getindex(count_t c) const
{
	if (c<5) return  sorted_[c];
	else return -1;
};

name_t site_t::getname(count_t c) const
{
	if (c<5){
		if (getcount(c)==0 ) return '*';
		return  profile::names_[sorted_[c]];
	}
	return '*';
};

name_t site_t::getname_gt(count_t c) const
{
	if (c<5){
		if (getcount(c)==getcount(c+1) || getcount(c)==0 ) return '*';
		return  profile::names_[sorted_[c]];
	}
	return '*';
};

std::string profile::getids(void)
{
	std::string str=decodeid0(site_.id0)+'\t'+decodeid1(site_.id1);
	for (unsigned int x=0; x<site_.extraid.size(); ++x){
		str+='\t'+decodeextraid(site_.extraid[x], x);
	};
	return str;
};

std::string profile::getids(const site_t &site)
{
	std::string str=decodeid0(site.id0)+'\t'+decodeid1(site.id1);
	for (unsigned int x=0; x<site.extraid.size(); ++x){
		str+='\t'+decodeextraid(site.extraid[x], x);
	};
	return str;
};
	
count_t profile::size(void) const
{
	return samples_;
};

count_t site_t::maskedcount(void) const
{
	count_t count=0;
	for(std::vector<quartet>::const_iterator it = sample.begin(); it != sample.end(); ++it) if (it->masked) count++;
	return count;
};

void profile::maskall(void){
	site_.maskall();
};
void site_t::maskall(void){
	for (unsigned int s=0; s<samples_;++s){
		sample[s].masked=true;
	};
};

void site_t::unmaskall(void){
	for (unsigned int s=0; s<samples_;++s){
		sample[s].masked=false;
	}
}

void profile::unmask(count_t a){
	if(a<samples_) site_.unmask(a);
}

void site_t::unmask(count_t a){
	if(a<samples_) sample[a].masked=false;
}

void profile::unmask(quartet_t *q){
	site_.unmask(q);
}

void site_t::unmask(quartet_t *q){
	q->masked=false;
}

void profile::mask(quartet_t *q){
	site_.mask(q);
}

void site_t::mask(quartet_t *q){
	q->masked=true;
}

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

void profile::setbasecount(const count_t &x, const count_t &y, const count_t &z){
	site_.sample[x].base[y]=z;
}
count_t count(const quartet_t c){
	return c.base[0]+c.base[1]+c.base[2]+c.base[3];
}

const uint64_t profile::getlinenumber(void) const{
	return size_;
}

const count_t profile::getid0(void) const{return site_.id0;}
const uint64_t profile::getid1(void) const{return site_.id1;}
void profile::setid0(const count_t &id0) {site_.id0=id0;}
void profile::setid1(const uint64_t &id1) {site_.id1=id1;}

const count_t profile::getextraid(const count_t &x) const {
	if (site_.extraid.size()>x) return site_.extraid[x];
	else return -1;
}

void profile::setextraid(const count_t &eid, const count_t &x) {
	if (site_.extraid.size()>x) {
		site_.extraid[x]=eid;
	} else {
		while (site_.extraid.size()<=x ) site_.extraid.push_back(-1);
		site_.extraid[x]=eid;
	}
}
	
void static merge(std::list <profile *>){
	
}
