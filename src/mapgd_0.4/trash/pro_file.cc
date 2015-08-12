#include "pro_file.h"

void profile::set_samples(const size_t &samples)
{
	sample_gof_=std::vector<real_t>(samples, 0);
	locus_->resize(samples);
}			

void profile::set_columns(const size_t &x) {
	columns_=x;
}

const std::string profile::get_sample_name(const size_t &a) const{
	return locus_->get_name(a);
}

int profile::set_sample_name(const std::string &str){
	for (size_t a; a<samples_; ++a) locus_->set_name(a, str+std::to_string(a) );
	return 0;
}

int profile::set_sample_name(const size_t &a, const std::string &str){
	return locus_->set_name(a, str);
}

int profile::read(row &this_row)
{
	if (mpileup_) return readm(this_row);
	else return readt(this_row);
}

int profile::readm(row &this_row) 
{
	std::string line;
	std::vector <std::string> column, quartet;

	if( std::getline(*in, line) ){
		column = split(line, delim_column);
		if(column.size() < 5){
			std::cerr << "WARNING:" << __FILE__ << ":" << __LINE__ << ": warning: skipping line ? with only " << column.size() << " fields." << std::endl;
			return 0;
		}

		std::vector <quartet_t>::iterator it=locus_->begin();
		std::vector <quartet_t>::iterator end=locus_->end();
		std::vector <std::string>::iterator column_it=column.begin();	

		id0_=index_.get_id0(column[0]);
		id1_=atoi(column[1].c_str());
		*ref_=(gt_t)(column[2][0]);

		column_it+=4;
		while (it!=end){
    			memset(it->base, 0, sizeof(count_t)*5 );
			scan(*column_it, *it);
			column_it+=3;
			++it;
		}
		return 0;
		
	} else {
		return EOF;
	}
	return UNEXPECTED;
}


///a nice function written by Berhard that scans the a line for quartets.
void inline profile::scan(const std::string &str, quartet_t &q)
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
			q.base[*ref_]++;
		} else {
			q.base[(gt_t)(*it)]++;
		}
		it++;
	}
}

///readt is a member of profile, which serves to convert between the three evil .pro formats, to the one nice clean map format.
int profile::readt(row &this_row){
	//The order of nucleotides read from a quartet file. Replace this with a memcopy of a constant.
	std::string line;

	std::vector <std::string> column, quartet;

	if (in==NULL){	//check to make sure that the file is open.
		std::cerr << __FILE__ << ":" << __LINE__ << ": error: attempted read line when no istream was open." << std::endl;
		exit(UNEXPECTED);
	};

	this_row.get("LOCUS", locus_);
	this_row.get("REF", ref_);
	this_row.get("ROWID", rowid_);

	std::vector<quartet_t>::iterator it=locus_->begin();

	if (std::getline(*in, line) ){	//read a line..
		switch (columns_){
		//READ A PRO FILE
			case 5:
			/* check for start of new scaffold*/
				if (line[0]!='>'){
					column=split(line, delim_column);
					if ( (column.size()-1) != samples_){
						std::cerr << __FILE__ << ":" << __LINE__ << ": error: could not parse line : \"" << line << "\"" << std::endl;
						std::cerr << column.size() << " fields found." << std::endl;
						std::cerr << "delimiter : \'" << delim_column << "\'" << std::endl;
						std::cerr << "columns : \'" << columns_ << "\'" << std::endl;
						exit(0);
						return EOF;
					}
					id1_=atoi(column[0].c_str() );
					for (unsigned int x=0; x< samples_; ++x){
						quartet=split(column[x+1], delim_quartet );
						it->base[0]=atoi(quartet[0].c_str() );
						it->base[1]=atoi(quartet[1].c_str() );
						it->base[2]=atoi(quartet[2].c_str() );
						it->base[3]=atoi(quartet[3].c_str() );
						++it;
					};
				} else {
					if(index_.open() ){
						index_.set_size(id0_, id1_);
						line.erase(line.begin());
						++id0_;
						index_.set_string(id0_, line);
					};
					read(this_row);
				}
				break;
			case 6:
				column=split(line, delim_column);
				if ( (column.size()-2) != samples_){
					std::cerr << __FILE__ << ":" << __LINE__ << ": error: could not parse line : \"" << line << "\"" << std::endl;
					std::cerr << column.size() << " fields found." << std::endl;
					std::cerr << "delimiter : \'" << delim_column << "\'" << std::endl;
					std::cerr << "columns : \'" << columns_ << "\'" << std::endl;
					exit(0);
					return EOF;
				};
	
				id1_=atoi(column[1].c_str() );
				id0_=index_.get_id0(column[0]);
				if(id0_!=last_id0_){
					index_.set_size(id0_, id1_);
					index_.set_string(id0_, column[0]);
					last_id0_=id0_;
				}
				for (size_t x=0; x < samples_; ++x){
					quartet=split(column[x+2], delim_quartet);
					it->base[0]=atoi(quartet[0].c_str() );
					it->base[1]=atoi(quartet[1].c_str() );
					it->base[2]=atoi(quartet[2].c_str() );
					it->base[3]=atoi(quartet[3].c_str() );
					++it;
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
				id1_=atoi(column[1].c_str() );
				id0_=index_.get_id0(column[0]);
				if(id0_!=last_id0_){
					index_.set_size(id0_, id1_);
					index_.set_string(id0_, column[0]);
					last_id0_=id0_;
				}
				*ref_=(gt_t)(column[2][0]);
				for (size_t x=0; x<samples_; ++x){
					quartet=split(column[x+3], delim_quartet);
					it->base[0]=atoi(quartet[0].c_str() );
					it->base[1]=atoi(quartet[1].c_str() );
					it->base[2]=atoi(quartet[2].c_str() );
					it->base[3]=atoi(quartet[3].c_str() );
					++it;
				};
				break;
		}
		return 0;
	}
	return EOF;
}

int profile::write(row &this_row){
	size_+=1;
	return writet(this_row);
}

int profile::writet(row &this_row){
	if (out==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempted write when no outstream was open." << std::endl;
		exit(UNEXPECTED);
	};

	/*these should really just be called in the constructor*/

	this_row.get("LOCUS", locus_);
	this_row.get("REF", ref_);
	this_row.get("ROWID", rowid_);

	id0_=index_.get_id0(*rowid_);
	id1_=index_.get_id1(*rowid_);

	std::vector<quartet_t>::iterator it=locus_->begin();
	std::vector<quartet_t>::iterator end=locus_->end();

	switch (columns_){
	case 5:
		if (id0_!=last_id0_){
			*out << '>' << index_.get_string(id0_) << std::endl;
		};
		*out << id1_;
		while(it!=end){
			if (it->masked){
				*out << delim_column;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0;
			}
			else {
				*out << delim_column;
				*out << it->base[0] << delim_quartet;
				*out << it->base[1] << delim_quartet;
				*out << it->base[2] << delim_quartet;
				*out << it->base[3];
			}
			++it;
		};
		*out << std::endl;
	break;
	case 6:
		*out << index_.get_string(id0_) << delim_column << id1_;
		*out << std::endl;
		while(it!=end){
			if (it->masked){
				*out << delim_column;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0;
			}
			else {
				*out << delim_column;
				*out << it->base[0] << delim_quartet;
				*out << it->base[1] << delim_quartet;
				*out << it->base[2] << delim_quartet;
				*out << it->base[3];
			}
			++it;
		}
	break;
	case 7: 
		*out << index_.get_string(id0_) << delim_column << id1_ << delim_column << (char)(*ref_);
		while(it!=end){
			if (it->masked){
				*out << delim_column;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0 << delim_quartet;
				*out << 0;
			}
			else {
				*out << delim_column;
				*out << it->base[0] << delim_quartet;
				*out << it->base[1] << delim_quartet;
				*out << it->base[2] << delim_quartet;
				*out << it->base[3];
			}
			++it;
		}
		*out << std::endl;	
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

	open_=false;
	in=NULL;
	out=NULL;

	if ( mode & std::fstream::in ){
		file.open(filename, std::ifstream::in);
		if (!file.is_open() ){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open " << filename << " for reading (1)." << std::endl;				
			exit(0);
		};
		in=&file;
		if (readheader()==BADHEADER){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot read header on " << filename << " (1). " << std::endl;
			std::cerr << "Vesions of mapgd >=0.3 require headers on .pro files" << std::endl;
			exit(0);
		};
		read_=true;
	} else if ( mode & std::fstream::out ){
		file.open(filename, std::ofstream::out);
		if (!file.is_open() ){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open " << filename << " for writing." << std::endl;				
			exit(0);
		};
		out=&file;
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
}

/** @brief closes a .pro file and unsets members.
  * @no return value
**/

void profile::close(void){
	if (file.is_open() ) file.close();
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
	
};

bool profile::is_open(void) const{
	return open_;
}

void profile::sort(void){
	locus_->sort();
}

size_t profile::size(void) const
{
	return samples_;
}

void profile::set_delim_quartet(const char &del)  
{
	delim_quartet=del;
}	

void profile::set_delim_column(const char &del)  
{
	delim_column=del;
}	

void profile::set_base_count(const size_t &x, const gt_t &y, const count_t &z)
{
	locus_->get_quartet(x).base[y]=z;
}

/*
gt_t profile::get_ref_(void) const 
{
	return ref_;
}

void profile::set_ref_(const char &eid)
{
	ref_=gt_t(eid);
}
*/
