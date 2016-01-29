#include "file-index.h"

const std::string File_index::file_name=".idx";
const std::string File_index::table_name="SCAFFOLDS";

File_index::File_index()
{
	cumulative_size_.push_back(0);
	id0_t last_id0_=-1;						
	std::string last_id0_str_="";
	open_=false;
}

File_index & File_index::operator =(const File_index& arg){
        id0_str_=arg.id0_str_;
	cumulative_size_=arg.cumulative_size_;
        id0_=arg.id0_;
        size_=arg.size_;
        last_id0_=arg.last_id0_;
        last_id0_str_=arg.last_id0_str_;
	return *this;
}


id1_t File_index::get_rowid (const id0_t &id0, const id1_t &id1) const{
	if (id1<=get_size(id0) ) return get_cumulative_size(id0)+id1;
	std::cerr << __FILE__ << ":" << __LINE__ << ": attempt to access positions outside of scaffold.\n";
	exit(0);
	return get_cumulative_size(id0)+id1;
}	

id1_t File_index::get_rowid (const std::string &id0_str, const id1_t &id1) const{
	return get_rowid(get_id0(id0_str), id1);
}


std::string File_index::get_string(const id0_t &id0) const{
	return id0_[id0];
}

id0_t File_index::get_id0 (const std::string &id0_str) const{
	std::map <std::string,id0_t>::const_iterator it;
	it=id0_str_.find(id0_str);
	if ( it!=id0_str_.end() ) return it->second;
	std::cerr << __FILE__ << ":" << __LINE__ << ":" << " no scaffold with name=" << id0_str << " exists.\n";
	exit(0);
	return ID0_MAX;
}

bool File_index::is_open (void) const{
	return open_;
};

id1_t File_index::get_size (const id0_t &id0) const{
	if ( id0<size_.size() ) return size_[id0];	
	std::cerr << __FILE__ << ":" << __LINE__ << ": scaffold does not exist.\n";
	exit(0);
	return ID0_MAX;
};

id1_t File_index::get_cumulative_size (const id0_t &id0) const{
	if ( id0<cumulative_size_.size() ) return cumulative_size_[id0];
	std::cerr << __FILE__ << ":" << __LINE__ << ":" << " no scaffold with id0=" << id0 << " exists.\n";
	exit(0);
	return ID0_MAX;
};

id1_t File_index::get_cumulative_size (const std::string &id0_str) const{
	if ( get_id0(id0_str)<cumulative_size_.size() ) return cumulative_size_[get_id0(id0_str)];
	std::cerr << __FILE__ << ":" << __LINE__ << ": scaffold does not exist.\n";
	return ID0_MAX;
};


std::istream& File_index::from_sam_header(std::istream &in)
{
	id1_t length;
	std::string line, val, tag;
	if (&in==NULL) {
		std::cerr << __FILE__ << ":" << __LINE__ << ": attempt to read from istream when no istream is open.\n";
		return in;
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
					cumulative_size_.push_back(cumulative_size_.back()+length);
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
	return in;
}

int File_index::write_index(std::ostream &out)
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

std::vector <id1_t> File_index::get_sizes (void)
{
	return size_;				
}	

std::ostream& operator << (std::ostream& out, const File_index& f)
{
	char delim_='\t';
	if (f.id0_.size()>0) out << f.id0_[0] << delim_ << f.size_[0];
	for (size_t x=1; x<f.id0_.size(); ++x){
		out << std::endl << f.id0_[x] << delim_;
		out << f.size_[x];
	}
	return out;
}

std::istream& operator >> (std::istream& in, File_index& f)
{
	std::string scaffold_name; 
	id1_t length;
//	char endl;
	in >> scaffold_name;
	while (scaffold_name!="@END_TABLE"){
		in >> length;
		f.size_.push_back(length);
		f.cumulative_size_.push_back(f.cumulative_size_.back()+length);
		f.id0_str_[scaffold_name]=f.id0_.size();
		f.id0_.push_back(scaffold_name);
		in >> scaffold_name;
	}
	in.get();
	return in;
}

std::string File_index::header(void) const
{
	return "@NAME          \tLENGTH\n";
}
	
size_t File_index::size(void) const
{
	return 2;
}
