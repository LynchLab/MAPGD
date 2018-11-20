#include "map_file.h"

//std::map <std::string, Map_file*(*)(void) > m_data_ctor;

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
	compressed_=false;
	indexed_=false;
}

bool 
Base_file::indexed(void) const
{
	return indexed_;
}

bool Base_file::check_compressed(void)
{
	if(open_) {
		return compressed_;
	} else  {
		uint16_t bytes;
		char c;
		in_->read((char *)(&bytes), 2);
		compressed_=(bytes==0x1f8b || bytes==0x8b1f);
		//std::cerr << "compressed?" << compressed_ << std::endl;
		return compressed_;
	}
}

const char & Base_file::get_delim(void) const
{
	return delim_column_;
}

bool Base_file::check_concatenated(void)
{
	if(open_) {
		return concatenated_;
	} else  {
		std::string line;
		std::vector <std::string> columns;
		std::getline(*in_, line);
		columns=split(line, '\t');
		concatenated_=std::find(columns.begin(), columns.end(), "CONCATENATED")!=columns.end();
		return concatenated_;
	}
}

void Base_file::open_no_extention(const char* filename, const std::ios_base::openmode &mode)
{
	if ( filename_.size()==0 ) filename_=std::string(filename);
#ifdef DEBUG
	std::cerr << __LINE__ << "opening no extension "<< filename <<"\n";
#endif 
	if ( open_ ){
		fprintf(stderr, gettext("mapgd:%s:%d: %s is already open.\n"), __FILE__, __LINE__, typeid(this).name() );
		exit(0);
	}
	if ( mode & std::ios::in ){
		file_.open( filename, std::ios::in);
		if ( !file_.is_open() )
		{
			fprintf(stderr, gettext("mapgd:%s:%d: cannot open %s for reading.\n"), __FILE__, __LINE__, filename);
			fprintf(stderr, gettext("Error: %s.\n"), strerror(errno) );
			exit(0);
		}
		buffer_.open(&in_, &file_);
		//in_=&file_;
		//std::cerr << "checking for compression...\n";
		buffer_.buffer_on();
		buffer_.reread_off();
		compressed_=check_compressed();
		if (compressed_)
		{
			//std::cerr << "Compressed!\n";
			file_.close();
			gzin_.open(filename, READ);
			if ( !gzin_.good() ){
				fprintf(stderr, gettext("mapgd:%s:%d: cannot open %s for reading.\n"), __FILE__, __LINE__, filename);
				fprintf(stderr, gettext("Error: %s.\n"), strerror(errno) );
				exit(0);
			};
			buffer_.open(&in_, &gzin_);
			buffer_.clear_read();
			buffer_.buffer_on();
			buffer_.reread_off();
		}
//		std::cerr << "checking for concatenation...\n";
		concatenated_=check_concatenated();
//		if(concatenated_) std::cerr << "Concatenated!\n";
		buffer_.buffer_off();
		buffer_.reread_on();
		read_=true;
		openmode_=mode;
	} else if ( mode & std::ios::out ){
		if (compressed_) 
		{
			gzout_.open( filename, std::ios::out);
			if (!gzout_.good() ){
				fprintf(stderr, gettext("mapgd:%s:%d: cannot open %s for writing.\n"), __FILE__, __LINE__, filename);
				fprintf(stderr, gettext("Error: %s.\n"), strerror(errno) );
				exit(0);
			};
			out_=&gzout_;
		} else {
			file_.open( filename, std::ios::out);
			if (!file_.is_open() ){
				fprintf(stderr, gettext("mapgd:%s:%d: cannot open %s for writing.\n"), __FILE__, __LINE__, filename);
				fprintf(stderr, gettext("Error: %s.\n"), strerror(errno) );
				exit(0);
			};
			out_=&file_;
		}
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
	else 
	{
		fprintf(stderr, gettext("mapgd:%s:%d: Unrecognized std::ios_base::openmode mode; cannot open file.\n"), __FILE__, __LINE__ );
	}
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

//         bool indexed(void) const;
//    bool concatenated(void) const;
bool 
Base_file::concatenated(void) 
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
//		std::cerr << "Base_file read  <" << data->table_name << ">" << char(in_->peek()) << std::endl;
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
//			std::cerr << "Base_file read (file_index) <" << data->table_name << ">" << char(in_->peek()) << std::endl;
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
Data *
Base_file::read_header(void)
{
	if (!read_ || !open_) 
	{
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to read from unopened file. Exiting.\n"), __FILE__, __LINE__ );
		exit(NOTOPEN);
	}

	std::string line;
	std::vector <std::string> columns, pair;
	std::getline(*in_, line);
	columns=split(line, '\t');
	if (columns.size()>2){
		pair=split(columns[0], ':');
		if (pair.size()==2) {
			if (pair[0]!="@NAME"){
				fprintf(stderr, gettext("mapgd:%s:%d: The first fields is not \"@NAME\". File cannot open properly. Exiting.\n"), __FILE__, __LINE__ );
				exit(BADHDR);
			}
			binary_=std::find(columns.begin(), columns.end(), "FORMAT:BINARY")!=columns.end();
#if (DEBUG)
			std::cerr << (columns[2]!="FORMAT:TEXT") << std::endl;
#endif
			concatenated_=std::find(columns.begin(), columns.end(), "CONCATENATED")!=columns.end();
			indexed_=std::find(columns.begin(), columns.end(), "INDEXED")!=columns.end();
			std::getline(*in_, line);
			columns=split(line, '\t');
			Data *data=Data::new_from_str(pair[1], columns);
			table_open_=true;
#if (DEBUG)
			std::cerr << "Successful, returning " << data->table_name << std::endl;
#endif
			return data;
		} else {
			fprintf(stderr, gettext("mapgd:%s:%d: The first fields is not a \"@A:B\" pair. File cannot open properly. Exiting.\n"), __FILE__, __LINE__ );
			exit(BADHDR);
		}
	}
	close();
	//This seems dangerous to me. The std::in doesn't seem to be reporting that it is eof when it is eof, thus, I am infering that it is 
	//eof from its failure to return a line. I really may need to go over to poll/status  
	if (line=="") return NULL;
	fprintf(stderr, gettext("mapgd:%s:%d: The following line was encountered:\n"), __FILE__, __LINE__ );
	std::cerr << line << std::endl;
	fprintf(stderr, gettext("mapgd:%s:%d: Header not terminated. File cannot open properly. Exiting.\n"), __FILE__, __LINE__ );
	exit(BADHDR);
}

bool 
Base_file::eof(void)
{
	if (open_){
		if (read_) return in_->eof();
		else if (write_) return out_->eof();
		return true;
	}
	else return true;
}

bool 
Base_file::binary(void) const
{
	return binary_;
}

void 
Base_file::write_text(File_index &file_index, const Indexed_data *data)
{
	*out_ << file_index.get_string(file_index.get_id0(data->get_abs_pos()) ) << '\t' << file_index.get_id1(data->get_abs_pos() ) << '\t' << *data << std::endl;
}

void 
Base_file::write_text(const Data *data)
{
	*out_ << *data << std::endl;
}

void 
Base_file::write_binary(const Data *data)
{
	data->write_binary(*out_);
}

void 
Base_file::write_binary(const Indexed_data *data)
{
	data->write_pos(*out_);
	data->write_binary(*out_);
}

Base_file & 
Base_file::write(File_index &file_index, const Indexed_data *data)
{
	if (open_ && write_)
	{
		if (binary_) write_binary(data);
		else write_text(file_index, data);
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to write to unopened file. Exiting.\n"), __FILE__, __LINE__ );
		exit(NOTOPEN);
	}
	return *this;
}

Base_file & 
Base_file::write(const Data *data)
{
	if (open_ && write_)
	{
		if (binary_) write_binary(data);
		else write_text(data);
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to write to unopened file. Exiting.\n"), __FILE__, __LINE__ );
		exit(NOTOPEN);
	}
	return *this;
}

void
Base_file::write_header(const Data *data)
{
	if (open_ && write_)
	{
		*out_ << "@NAME:" << data->get_table_name() << "\tVERSION:" << VERSION;
		if (try_binary_ && data->get_binary() )
		{
			*out_ << "\tFORMAT:BINARY";
			binary_=true;
		} else {
			*out_ << "\tFORMAT:TEXT";
			binary_=false;
		}
		if (concatenated_) *out_ << "\tCONCATENATED";
			*out_ << std::endl;
		*out_ << data->header();
		table_open_=true;
	} else {
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to write to unopened file. Exiting.\n"), __FILE__, __LINE__ );
		exit(NOTOPEN);
	}
}

void
Base_file::write_header(const File_index &file_index, const Data *data)
{
	if (!write_ || !open_) 
	{
		fprintf(stderr, gettext("mapgd:%s:%d: Attempt to write to unopened file. Exiting.\n"), __FILE__, __LINE__ );
		exit(NOTOPEN);
	}

	Flat_file <File_index> index;

#ifdef DEBUG
	std::cerr << "Writing header of " << filename_ << std::endl;
	std::cerr << "aka: " << this->filename() << std::endl;
#endif
	index.open_from(*this);

	if (!index.is_open() ) 
	{
		fprintf(stderr, gettext("mapgd:%s:%d: Could not write to file. Exiting.\n"), __FILE__, __LINE__ );
		exit(NOTOPEN);
	}
	std::cerr << "writing header\n";
	index.write_header(file_index);
	std::cerr << "writing file\n";
	index.write(file_index);
	index.close_table();
	*out_ << "@NAME:" << data->get_table_name() << "\tVERSION:" << VERSION;
	if (try_binary_ && data->get_binary() ){
		*out_ << "\tFORMAT:BINARY";
		binary_=true;
	} else {
		*out_ << "\tFORMAT:TEXT";
		binary_=false;
	}
	if (concatenated_) *out_ << "\tCONCATENATED";
	*out_ << "\tINDEXED\n";
	*out_ << data->header();
	table_open_=true;
}
