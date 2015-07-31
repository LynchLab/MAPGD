#include "file-index.h"

file_index::file_index()
{
	id0_t last_id0_=-1;						
	std::string last_id0_str_="";
	char delim_='\t';				
	char key_='@';			
}

int file_index::from_sam_header(std::istream &in)
{
	id1_t length;
	std::string line, sn, ln;
	if (&in==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempt to read from istream when no istream is open.\n";
		return -1;
	};
	while ( std::getline(in, line) ) {
		std::vector <std::string> column=split(line, '\t');
		if ( column[0]=="@SQ" ) {
			if (column.size()==3) {
				sn=split_first(column[1], ':')[1];
				ln=split(column[2], ':')[1];
				last_id0_str_=sn;
				length=atoi(ln.c_str());
				last_id0_=id0_.size();
				id0_str_[last_id0_str_]=last_id0_;
				id0_.push_back(last_id0_str_);
				size_.push_back(length);
			} else {
				std::cerr << __FILE__ << ":" << __LINE__ << ": incorrect number of columns in samfile header.\n";
			}; 
		};
	}
	return 0;
}

int file_index::write_index(std::ostream &out)
{
	id1_t length;
	std::string line, sn, ln;
	if (&out==NULL){
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempt to write to ostream when no ostream is open.\n";
		return -1;
	}
	for (size_t x; x<id0_.size(); ++x){
		out << "@SQ" << delim_;
		out << "SN:" << id0_[x] << delim_; 
		out << "LN:" << size_[x] << std::endl; 
	}
	return 0;
}

id0_t file_index::encodeid0(const std::string &id)
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
