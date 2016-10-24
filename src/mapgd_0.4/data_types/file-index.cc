#include "file-index.h"

const std::string File_index::file_name=".idx";
const std::string File_index::table_name="SCAFFOLDS";
const Registration File_index::registered=Registration(File_index::table_name, File_index::create);

File_index::File_index()
{
	cumulative_size_.push_back(0);
	last_id0_=-1;						
	std::string last_id0_str_="";
	open_=false;
}

File_index & 
File_index::operator =(const File_index& arg)
{
        id0_str_=arg.id0_str_;
	cumulative_size_=arg.cumulative_size_;
        id0_=arg.id0_;
        size_=arg.size_;
        last_id0_=arg.last_id0_;
        last_id0_str_=arg.last_id0_str_;
	return *this;
}


id1_t 
File_index::get_abs_pos (const id0_t &id0, const id1_t &id1) const
{
	if (id1<=get_size(id0) ) return get_cumulative_size(id0)+id1;
	std::cerr << __FILE__ << ":" << __LINE__ << ": attempt to access positions outside of scaffold.\n";
	exit(0);
	return get_cumulative_size(id0)+id1;
}	

id1_t 
File_index::get_abs_pos (const std::string &id0_str, const id1_t &id1) const
{
	return get_abs_pos(get_id0(id0_str), id1);
}

std::string 
File_index::get_string(const id0_t &id0) const
{
	return id0_[id0];
}

id0_t 
File_index::get_id0 (const std::string &id0_str) const
{
	std::map <std::string,id0_t>::const_iterator it;
	it=id0_str_.find(id0_str);
	if ( it!=id0_str_.end() ) return it->second;
	std::cerr << __FILE__ << ":" << __LINE__ << ":" << " no scaffold with name=" << id0_str << " exists.\n";
	exit(0);
	return ID0_MAX;
}

id1_t
File_index::get_id1(const id1_t &abs_pos) const
{
	id0_t id0=get_id0(abs_pos);
//	std::cerr << abs_pos << ":" << id0 << "::" << get_cumulative_size(id0) << std::endl;
	return abs_pos-get_cumulative_size(id0);
}

id0_t
File_index::get_id0(const id1_t &abs_pos) const
{
	id0_t id0=1;
	while(abs_pos>cumulative_size_[id0]) id0++;
	return id0-1;
}

bool 
File_index::is_open (void) const
{
	return open_;
}

id1_t 
File_index::get_size (const id0_t &id0) const
{
	if ( id0<size_.size() ) return size_[id0];	
	std::cerr << __FILE__ << ":" << __LINE__ << ": scaffold does not exist.\n";
	exit(0);
	return ID0_MAX;
}

id1_t 
File_index::get_cumulative_size (const id0_t &id0) const
{
	if ( id0<cumulative_size_.size() ) return cumulative_size_[id0];
	std::cerr << __FILE__ << ":" << __LINE__ << ":" << " no scaffold with id0=" << id0 << " exists.\n";
	exit(0);
	return ID0_MAX;
}

id1_t 
File_index::get_cumulative_size (const std::string &id0_str) const
{
	if ( get_id0(id0_str)<cumulative_size_.size() ) return cumulative_size_[get_id0(id0_str)];
	std::cerr << __FILE__ << ":" << __LINE__ << ": scaffold does not exist.\n";
	return ID0_MAX;
}

id1_t
File_index::get_reference_size(void) const
{
	return cumulative_size_.back();
}


std::istream & 
File_index::from_sam_header(std::istream &in)
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

int 
File_index::write_index(std::ostream &out)
{
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

const std::vector <id1_t> 
File_index::get_sizes (void) const
{
	return size_;				
}

const std::vector<std::string>
File_index::get_names (void) const
{
	return id0_;				
}	

void
File_index::write(std::ostream& out) const
{
	char delim_='\t';
	if (id0_.size()>0) out << id0_[0] << delim_ << size_[0];
	for (size_t x=1; x<id0_.size(); ++x){
		out << std::endl << id0_[x] << delim_;
		out << size_[x];
	}
}

void 
File_index::read (std::istream& in)
{
	std::string scaffold_name; 
	id1_t length;
	std::string line;
	if (in.good() ){
		while ( in.peek()!='@' && in.good() ) {
			std::getline(in, line);
			std::stringstream stream_line(line);
			stream_line >> scaffold_name;
			stream_line >> length;
			size_.push_back(length);
			cumulative_size_.push_back(cumulative_size_.back()+length);
			id0_str_[scaffold_name]=id0_.size();
			id0_.push_back(scaffold_name);
		} 
	}
	
}

std::string 
File_index::header(void) const
{
	return "@NAME          \tLENGTH\n";
}
	
size_t 
File_index::size(void) const
{
	return 2;
}

const std::string 
File_index::sql_header(void) const {
	return "(SCFNAME varchar(255), START int, STOP int)";
}

const std::string 
File_index::sql_column_names(void) const {
	return "(SCFNAME, STOP, START)";
}

const std::string 
File_index::sql_values(void) const {
        char return_buffer[SQL_LINE_SIZE]={0};
	char *write_ptr=return_buffer;
	std::vector<std::string>::const_iterator id0_it=id0_.cbegin(), end=id0_.cend();
	std::vector<id1_t>::const_iterator size_it=cumulative_size_.cbegin();

	write_ptr+=snprintf(return_buffer, SQL_LINE_SIZE, "('%s','%d', '%d')", 
       		(id0_it++)->c_str(),
		*(++size_it),
		*size_it );
	while (id0_it!=end) {
		write_ptr+=snprintf(write_ptr, SQL_LINE_SIZE, ", ('%s','%d', '%d')", 
       		(id0_it++)->c_str(),
		*(++size_it),
		*size_it);
        }
	return std::string(return_buffer);
}

void
File_index::sql_read (std::istream &in) 
{	
	std::string scaffold_name;
	id1_t start, stop, length;
	in >> scaffold_name;
	in >> start;
	in >> stop;
	length=stop-start;
	size_.push_back(length);
	cumulative_size_.push_back(cumulative_size_.back()+length);
	id0_str_[scaffold_name]=id0_.size();
	id0_.push_back(scaffold_name);
	//just want to read till the new line, not doing much else.
	std::getline(in, scaffold_name);
}
