#include "map-file.h"

Base_file::Base_file(void)
{
	open_=false;
	table_open_=false;
	read_=false;
	write_=false;
	delim_column_='\t';
	binary_=false;
	try_binary_=false;
	filename_="";
	concatenated_=false;
	indexed_=false;
}

bool 
Base_file::indexed(void)
{
	return indexed_;
}

bool Base_file::check_concatenated(const char* filename)
{
	std::cerr << "opening\n";
	if(open_) {
		return concatenated_;
	} else  {
		if ( filename_.size()==0 ) filename_=std::string(filename);
		file_.open( filename, std::ios::in);
		if ( !file_.is_open() ){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open " << filename << " for reading." << std::endl;
			std::cerr << "Error: " << strerror(errno) << std::endl;
			exit(0);
		}
		in_=&file_;
		std::string line;
		std::vector <std::string> columns;
		std::getline(*in_, line);
		columns=split(line, '\t');
		concatenated_=std::find(columns.begin(), columns.end(), "CONCATENATED")!=columns.end();
		file_.close();
		in_=NULL;
		return concatenated_;
	}
}

void Base_file::open_no_extention(const char* filename, const std::ios_base::openmode &mode)
{
	if ( filename_.size()==0 ) filename_=std::string(filename);
#ifdef DEBUG
	std::cerr << "opening "<< filename <<"\n";
#endif 
	if ( open_ ){
		std::cerr << __FILE__ << ":" << __LINE__ << ": " << typeid(this).name() << " is already open." << std::endl;
		exit(0);
	}
	if ( mode & std::ios::in ){
		file_.open( filename, std::ios::in);
		if ( !file_.is_open() ){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open " << filename << " for reading." << std::endl;
			std::cerr << "Error: " << strerror(errno) << std::endl;
			exit(0);
		};
		in_=&file_;
		read_=true;
		openmode_=mode;
	} else if ( mode & std::ios::out ){
		file_.open( filename, std::ios::out);
		if (!file_.is_open() ){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open " << filename << " for writing." << std::endl;
			std::cerr << "Error: " << strerror(errno) << std::endl;
			exit(0);
		};
		out_=&file_;
		write_=true;
		openmode_=mode;
	};
	if ( mode & std::ios::binary ) {
		try_binary_=true;
		binary_=false;
	}
	open_=true;
}

void Base_file::open(const std::ios_base::openmode &mode)
{
	if ( mode & std::ios::in ) open(&std::cin, mode);
	else if ( mode & std::ios::out ) open(&std::cout, mode);
}


void Base_file::open(std::iostream *s, const std::ios_base::openmode &mode)
{
	#ifdef POSIX
	if ( check_stream(s) ){
		std::cerr << __FILE__ << ":" << __LINE__ << ": iostream is empty. Exiting." << std::endl;
		exit(0);
	}
	#endif
	if ( open_ ){
		std::cerr << __FILE__ << ":" << __LINE__ << ": " << typeid(this).name() << " is already open." << std::endl;
		exit(0);
	}
	if ( mode & std::ios::in ) {
		in_=s;
		read_=true;
		openmode_=mode;
	}
	if ( mode & std::ios::out ) {
		out_=s;
		write_=true;
		openmode_=mode;
	}
	if ( mode & std::ios::binary ) {
		try_binary_=true;
		binary_=true;
	}
	open_=true;
	concatenated_=true;
}

void Base_file::open(std::istream *s, const std::ios_base::openmode &mode)
{
#ifdef POSIX
	if (check_stream(s) ) {
		std::cerr << __FILE__ << ":" << __LINE__ << ": iostream is empty. Exiting." << std::endl;
		exit(0);
	}
#endif
	if ( open_ ) {
		std::cerr << __FILE__ << ":" << __LINE__ << ": " << typeid(this).name() << " is already open." << std::endl;
		exit(0);
	}
	if ( mode & std::ios::in) {
		in_=s;
		read_=true;
		openmode_=mode;
	}
	if ( mode & std::ios::binary) {
		binary_=false;
		try_binary_=true;
	}
	open_=true;
	concatenated_=true;
}

void Base_file::open(std::ostream *s, const std::ios_base::openmode &mode)
{
	if ( open_ ){
		std::cerr << __FILE__ << ":" << __LINE__ << ": " << typeid(this).name() << " is already open." << std::endl;
		exit(0);
	}
	if ( mode & std::ios::out) {
		out_=s;
		write_=true;
		openmode_=mode;
	}
	if ( mode & std::ios::binary){
		binary_=false;
		try_binary_=true;
	}
	open_=true;
	concatenated_=true;
}

std::istream * Base_file::get_in(void)
{
	return in_;
}

std::ostream * Base_file::get_out(void)
{
	return out_;
}

const std::ios::openmode& Base_file::openmode(void)
{
	return openmode_;
}

const bool& Base_file::concatenated(void)
{
	return concatenated_;
}

void Base_file::seekp(const std::streampos &pos)
{
	out_->seekp(pos);
}

void Base_file::seekg(const std::streampos &pos)
{
	in_->seekg(pos);
}

std::streampos Base_file::tellp(void)
{
	return out_->tellp();
}

std::streampos Base_file::tellg(void)
{
	return in_->tellg();
}

const std::string& Base_file::filename(void)
{
	return filename_;
}

void Base_file::close_table(void){
#ifdef DEBUG
	std::cerr << "Closing table " << filename_ << std::endl;
#endif
	if (write_ && open_) *out_ << "@END_TABLE\n";
	table_open_=false;
}

void Base_file::close(void)
{
	if (table_open_){
		close_table();
	}
	if ( file_.is_open() ) {
		file_.close();
	}
	open_=false;
}

Base_file & 
Base_file::read(Data *data)
{
	if (!table_open_){
		std::cerr << __FILE__<< ":" <<__LINE__ << ": file not open. The methods Base_file::open() and Base_file::read_header(Data *) should be called.";
		return *this;
	}
	if (read_){
		if (binary_){
			in_->read( (char *) data, data->size() );
		}
		else {
			if (in_->peek()=='@') {
				std::string line;
				std::getline(*in_, line);
				if (line!="@END_TABLE") {
					std::cerr << __FILE__ << ":" << __LINE__ << ": file not closed correctly, exiting.\n";
					exit(0);
				}
				if (!concatenated_) {
					this->close();
				} else {
					this->close_table();
				}
			} else *in_ >> *data;
#ifdef SAFER
			if (in_->peek()=='\n') {
				std::cerr << __FILE__<< ":" <<__LINE__ << ": line not parsed correctly.";
				in_->get();
			}
#endif
		}
	} else {
		std::cerr << __FILE__<< ":" <<__LINE__ << ": file not open. The methods Base_file::open() and Base_file::read_header(Data *) should be called.";
	}
	return *this;
}

Base_file &
Base_file::read(File_index &index, Indexed_data *data)
{
	std::string scaffold;
	id1_t pos;
	if (!table_open_){
		std::cerr << __FILE__<< ":" <<__LINE__ << ": file not open. The methods Base_file::open() and Base_file::read_header(Data *) should be called.";
		return *this;
	}
	if (read_){
		if (binary_){
			in_->read( (char *) data, data->size() );
		}
		else {
			if (in_->peek()=='@') {
				std::string line;
				std::getline(*in_, line);
				if (line!="@END_TABLE") {
					std::cerr << __FILE__ << ":" << __LINE__ << ": file not closed correctly, exiting.\n";
					exit(0);
				}
				if (!concatenated_) {
					this->close();
				} else {
					this->close_table();
				}
			} else {
				*in_ >> scaffold;
				*in_ >> pos;
				*in_ >> *data;
				data->set_abs_pos(index.get_abs_pos(scaffold, pos) );
			}
#ifdef SAFER
			if (in_->peek()=='\n') {
				std::cerr << __FILE__<< ":" <<__LINE__ << ": line not parsed correctly.";
				in_->get();
			}
#endif
		}
	} else {
		std::cerr << __FILE__<< ":" <<__LINE__ << ": file not open. The methods Base_file::open() and Base_file::read_header(Data *) should be called.";
	}
	return *this;
}

bool Base_file::is_open(void) const
{
	return open_;
}			

bool Base_file::table_is_open(void) const
{
	return table_open_;
}			


//read_header needs to be split into two stages.
Data *Base_file::read_header(void)
{
	std::string line;
	std::vector <std::string> columns, pair;
	std::getline(*in_, line);
	columns=split(line, '\t');
	if (columns.size()>2){
		pair=split(columns[0], ':');
		if (pair.size()==2) {
			if (pair[0]!="@NAME"){
				std::cerr << __FILE__ << ":" << __LINE__ << ": file cannot open properly. Exiting.\n"; 
				exit(0);
			}
			binary_=std::find(columns.begin(), columns.end(), "FORMAT:BINARY")!=columns.end();
#if (DEBUG)
			std::cerr << (columns[2]!="FORMAT:TEXT") << std::endl;
#endif
			std::cerr << "Checking 3" << std::endl;
			concatenated_=std::find(columns.begin(), columns.end(), "CONCATENATED")!=columns.end();
			indexed_=std::find(columns.begin(), columns.end(), "INDEXED")!=columns.end();
			std::getline(*in_, line);
			columns=split(line, '\t');
			Data *data=Data::new_from_str(pair[1], columns);
			table_open_=true;
			return data;
		} else {
			std::cerr << __FILE__ << ":" << __LINE__ << ": file cannot open properly. Exiting.\n"; 
			exit(0);
		}
	}
	close();
	//This seems dangerous to me. The std::in doesn't seem to be reporting that it is eof when it is eof, thus, I am infering that it is 
	//eof from its failure to return a line. I really may need to go over to poll/status  
	if (line=="") return NULL;
	std::cerr << line << std::endl;
	std::cerr << __FILE__ << ":" << __LINE__ << ": file cannot open properly. Exiting.\n"; 
	exit(0);
}

bool Base_file::eof(void)
{
	if (open_){
		if (read_) return in_->eof();
		else if (write_) return out_->eof();
		return true;
	}
	else return true;
}

bool Base_file::binary(void) const
{
	return binary_;
}
