#include "file-index.h"

file_index::file_index()
{
	id0_t last_id0_=-1;						
	std::string last_id0_str_="";
	open_=false;
}

file_index & file_index::operator =(const file_index& arg){
        id0_str_=arg.id0_str_;
        id0_=arg.id0_;
        size_=arg.size_;
        last_id0_=arg.last_id0_;
        last_id0_str_=arg.last_id0_str_;
	return *this;
}


std::string file_index::get_string(const id0_t &id0) const{
	return id0_[id0];
}

id0_t file_index::get_id0 (const std::string &id0_str) const{
	std::map <std::string,id0_t>::const_iterator it;
	it=id0_str_.find(id0_str);
	if (it!=id0_str_.end() ) return it->second;
	else return NULL;
}

bool file_index::is_open (void) const{
	return open_;
};

id1_t file_index::get_size (const id0_t &id0) const{
	if (id0<size_.size() )return size_[id0];	
	else return 0;
};

int file_index::from_sam_header(std::istream &in)
{
	id1_t length;
	std::string line, val, tag;
	if (&in==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempt to read from istream when no istream is open.\n";
		return -1;
	};
	while ( std::getline(in, line) ) {
		std::vector <std::string> column=split(line, '\t');
		bool has_sn=false, has_ln=false;
		if ( column[0]=="@SQ" ) {
			for (size_t x=1; x<column.size(); ++x){
				tag=split(column[x], ':')[0];
				val=split(column[x], ':')[1];
				if(tag.compare("LN")==0){
					length=atoi(val.c_str());
					size_.push_back(length);
					has_ln=true;
				} else if(tag.compare("SN")==0) {
					last_id0_str_=val;
					last_id0_=id0_.size();
					id0_str_[last_id0_str_]=last_id0_;
					id0_.push_back(last_id0_str_);
					has_sn=true;
				}
			}
			if( not (has_sn && has_ln) ) std::cerr << __FILE__ << ":" << __LINE__ << ": formating error in samfile header.\n";
		};
	}
	open_=true;
	return 0;
}

int file_index::write_index(std::ostream &out)
{
	id1_t length;
	char delim_='\t';
	std::string line, sn, ln;
	if (&out==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempt to write to ostream when no ostream is open.\n";
		return -1;
	}
	for (size_t x=0; x<id0_.size(); ++x){
		out << "@SQ" << delim_;
		out << "SN:" << id0_[x] << delim_; 
		out << "LN:" << size_[x] << std::endl; 
	}
	return 0;
}

std::vector <id1_t> file_index::get_sizes (void)
{
	return size_;				
}	

