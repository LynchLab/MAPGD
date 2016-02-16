#include "map-file.h"

Base_file::Base_file(void)
{
	open_=false;
	table_open_=false;
	read_=false;
	write_=false;
	delim_column_='\t';
	binary_=false;
	filename_="";
	concatenated_=false;
	indexed_=false;
}

bool 
Base_file::indexed(void)
{
	return indexed_;
}


//This seriously needs to be cleaned up.
template <class T>
void Data_file<T>::open(const char* filename, const std::ios_base::openmode &mode)
{
	std::vector<std::string> line;
	filename_=std::string(filename);
	if (open_){
		std::cerr << __FILE__ << ":" << __LINE__ << ": " << typeid(this).name() << " is already open." << std::endl;
		exit(0);
	}
	if ( mode & std::ios::in ){
		line=split_last(filename, '.');
		if (line.back()!=T::file_name) {
#ifdef DEBUG
			std::cerr << line.back() << "!=" << T::file_name;
			std::cerr << "opening " << filename_ << " in concatenated mode\n";
#endif
			concatenated_=true;
			open_no_extention(filename, mode);
			return;
		} else {
			concatenated_=false;
			filename_=line[0];
		}
	} else {
		concatenated_=false;
	} 
#ifdef DEBUG
	std::cerr << "opening "<< filename_ <<" in split mode\n";
#endif 
	open_extention(filename_.c_str(), mode);
}

template <class T>
void Data_file<T>::open_extention(const char* filename, const std::ios_base::openmode &mode)
{
	if (filename_.size()==0) filename_=std::string(filename);
	std::string temp_filename=std::string(filename)+T::file_name;
#ifdef DEBUG
	std::cerr << "opening "<< temp_filename <<"\n";
#endif 
	open_no_extention(temp_filename.c_str(), mode);
}

void Base_file::open_no_extention(const char* filename, const std::ios_base::openmode &mode)
{
	if (filename_.size()==0) filename_=std::string(filename);
#ifdef DEBUG
	std::cerr << "opening "<< filename <<"\n";
#endif 
	if (open_){
		std::cerr << __FILE__ << ":" << __LINE__ << ": " << typeid(this).name() << " is already open." << std::endl;
		exit(0);
	}
	if ( mode & std::ios::in ){
		file_.open( filename, std::ios::in);
		if (!file_.is_open() ){
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
	if (mode & std::ios::binary) binary_=true;
	open_=true;
}

void Base_file::open(const std::ios_base::openmode &mode)
{
	if ( mode & std::ios::in ) open(&std::cin, mode);
	else if ( mode & std::ios::out) open(&std::cout, mode);
}


void Base_file::open(std::iostream *s, const std::ios_base::openmode &mode)
{
	#ifdef POSIX
	if (check_stream(s) ){
		std::cerr << __FILE__ << ":" << __LINE__ << ": iostream is empty. Exiting." << std::endl;
		exit(0);
	}
	#endif
	if ( open_ ){
		std::cerr << __FILE__ << ":" << __LINE__ << ": " << typeid(this).name() << " is already open." << std::endl;
		exit(0);
	}
	if ( mode & std::ios::in) {
		in_=s;
		read_=true;
		openmode_=mode;
	}
	if ( mode & std::ios::out) {
		out_=s;
		write_=true;
		openmode_=mode;
	}
	if ( mode & std::ios::binary) binary_=true;
	open_=true;
	concatenated_=true;
}

void Base_file::open(std::istream *s, const std::ios_base::openmode &mode)
{
#ifdef POSIX
	if (check_stream(s) ){
		std::cerr << __FILE__ << ":" << __LINE__ << ": iostream is empty. Exiting." << std::endl;
		exit(0);
	}
#endif
	if ( open_ ){
		std::cerr << __FILE__ << ":" << __LINE__ << ": " << typeid(this).name() << " is already open." << std::endl;
		exit(0);
	}
	if ( mode & std::ios::in) {
		in_=s;
		read_=true;
		openmode_=mode;
	}
	if ( mode & std::ios::binary) binary_=true;
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
	if ( mode & std::ios::binary) binary_=true;
	open_=true;
	concatenated_=true;
}

template <class T>
void Data_file<T>::open_from(Base_file &file)
{
	if (file.table_is_open() ) file.close_table();
	if (file.openmode() & std::ios::in){
//		if(file.filename().size()!=0) {
		if(file.concatenated() ) this->open(file.get_in(), file.openmode() );
		else this->open_extention(file.filename().c_str(), file.openmode() );
//		}
//		else this->open(std::ios::in);
	} else if (file.openmode() & std::ios::out) {
//		if(file.filename().size()!=0) {
		if(file.concatenated() ) this->open(file.get_out(), file.openmode() );
		else this->open_extention(file.filename().c_str(), file.openmode() );
//		}
//		else this->open(std::ios::out);
	}
	binary_=(file.openmode() & std::ios::binary);
	open_=true;
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

/* May replace open header
void Base_file::open_table(void)
{
	table_open_=true;
}*/

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

template <class T>
id1_t Indexed_file<T>::get_pos(const T &data) const 
{
	return data.get_abs_pos();
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

template <class T>
Data_file<T>& Data_file<T>::read(T &data)
{
	if (!table_open_ ){
		return *this;
	}
	if (read_){
		if (binary_) read_binary(data);
		else read_text(data);
	} else {
		std::cerr << __FILE__<< ":" <<__LINE__ << ": file not open for reading. The methods Flat_file<type>::open() and Flat_file<type>::read_header(<type>) should be called.";
	}
	return *this;
}

template <class T>
void Flat_file<T>::read_text(T &data)
{
	if (!in_->good() ) {
		std::cerr << "an error has occured durring reading.\n";
		exit(0);
	}
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
	}
	else *in_ >> data;
}

template <class T>
void Indexed_file<T>::read_text(T &data)
{
	//TODO Check for table_open instead?
	id1_t pos;
	std::string scaffold;
	if (!in_->good() ) {
		std::cerr << "an error has occured durring reading.\n";
		exit(0);
	}
	if (in_->peek()=='@') {
		std::string line;
		*in_ >> line;
		if (line!="@END_TABLE") {
			std::cerr << line << std::endl;
			std::cerr << __FILE__ << ":" << __LINE__ << ": file not closed correctly, exiting.\n";
			exit(0);
		} 
		if (!concatenated_) {
			this->close();
		} else {
			this->close_table();
		}
	} else {
#ifdef DEBUG
		std::cerr << (char)(in_->peek()) << std::endl;
#endif
		*in_ >> scaffold;
		*in_ >> pos;
		*in_ >> data;
		data.set_abs_pos(file_index_.get_abs_pos(scaffold, pos) );
	}
}

template <class T>
void Mpileup_file<T>::read_text(T &data)
{
	//TODO Check for table_open instead?
	std::string scaffold;
	id1_t pos;
	*in_ >> scaffold;
	*in_ >> pos;
	if (!in_->good() ) {
		this->close();
	} else {
		data.set_abs_pos(file_index_.get_abs_pos(scaffold, pos) );
		mpileup(*in_, data);
		if (in_->eof() ) this->close();
	}
}

template <class T>
void Data_file<T>::read_binary(T &data)
{
	in_->read( (char *) &data, data.size() );
}

template <class T>
void Flat_file<T>::write_text(const T &data)
{
	*out_ << data << std::endl;
}

template <class T>
void Indexed_file<T>::write_text(const T &data)
{
	*out_ << file_index_.get_string(file_index_.get_id0(data.get_abs_pos()) ) << '\t' << file_index_.get_id1(data.get_abs_pos() ) << '\t' << data << std::endl;
}


template <class T>
void Data_file<T>::write_binary(const T &data)
{
	out_->write( (const char *) &data, data.size() );
}

template <class T>
Data_file<T>& Data_file<T>::write(const T &data)
{
	if (binary_) write_binary(data);
	else write_text(data);
//	if (!out_->good() ) { std::cerr << __FILE__ << ":" << __LINE__ << ": unexpected error writing file. Exiting.\n"; exit(0);};
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

template <class T>
void Flat_file<T>::write_header(const T &data)
{
	if (write_) {
		*out_ << "@NAME:" << T::table_name << "\tVERSION:" << VERSION;
		if (binary_) *out_ << "\tFORMAT:BINARY";
		else *out_ << "\tFORMAT:TEXT";
		if (concatenated_) *out_ << "\tCONCATENATED";
		*out_ << std::endl;
		*out_ << data.header();
		table_open_=true;
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ << ": file not open for writing. Exiting.\n"; exit(0);
	}
}

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
			binary_=(columns[2]!="FORMAT:TEXT");
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

template <class T>
T Flat_file<T>::read_header(void)
{
	std::string line;
	std::vector <std::string> columns;
	std::getline(*in_, line);
	columns=split(line, '\t');
	if (columns.size()>2){
		if (columns[0]=="@NAME:"+T::table_name ){
			binary_=(columns[2]!="FORMAT:TEXT");
			std::getline(*in_, line);
			columns=split(line, '\t');
			T data(columns);
			table_open_=true;
			return data;
		}
		std::cerr << __FILE__ << ":" << __LINE__ << " attempted to open incorrect header.\n"; 
	}
	table_open_=false;
	std::cerr << __FILE__ << ":" << __LINE__ << " could not initilize " << typeid(T).name() <<"\n";
	T data;
	return data;
}

template <class T>
T Indexed_file<T>::read_header(void)
{
	Flat_file <File_index> index;
	index.open_from(*this);
	file_index_=index.read_header();
	while(index.table_is_open() ){
		index.read(file_index_);
	}
	
	std::string line;
	std::vector <std::string> columns;
#ifdef DEBUG
	std::cerr << in_->peek() << std::endl;
#endif
	/*while(in_->peek()!='@' && in_->good() ){
		std::getline(*in_, line);
		std::cerr << line << std::endl;
	}*/
	std::getline(*in_, line);
	
	columns=split(line, '\t');
	if (columns.size()>2){
		if (columns[0]=="@NAME:"+T::table_name){
			//binary_=(columns[2]!="FORMAT:TEXT");
			std::getline(*in_, line);
			columns=split(line, '\t');
			T data(columns);
			table_open_=true;
			return data;
		}
		std::cerr << __FILE__ << ":" << __LINE__ << " attempted to open incorrect header.\n"; 
		std::cerr << line << std::endl;
		std::cerr << T::table_name << std::endl;
	}
	std::cerr << __FILE__ << ":" << __LINE__ << " could not initilize " << typeid(T).name() << "\n";
	std::cerr << line << std::endl;
	T data;
	table_open_=false;
	return data;
}

template <class T>
T Mpileup_file<T>::read_header(void)
{
	std::string line;
	std::vector <std::string> columns;
	std::getline(*in_, line);
	columns=split(line, '\t');
	T data( (columns.size()-3.)/3. );	//TODO FIX IMEDIATELY!! (TOMORROW).
	in_->putback('\n');
	for(std::string::iterator c=line.end(); (c--)!=line.begin();) {
		in_->putback(*c);
	}
	table_open_=true;
	return data;
}

template <class T>
void Indexed_file<T>::write_header(const T &data)
{
	if (!write_) {
		std::cerr << __FILE__ << ":" << __LINE__ << " file not open for writing. Exiting \n";
		exit(0);
	}
	Flat_file <File_index> index;
#ifdef DEBUG
	std::cerr << "Writing header of " << filename_ << std::endl;
	std::cerr << "aka: " << this->filename() << std::endl;
#endif
	index.open_from(*this);
	if (!index.is_open() ) {
		std::cerr << __FILE__ << ":" << __LINE__ << " cannot  open for writing. Exiting \n";
		exit(0);
	}
	index.write_header(file_index_);
	index.write(file_index_);
	index.close_table();
	*out_ << "@NAME:" << T::table_name << "\tVERSION:" << VERSION;
	if (binary_) *out_ << "\tFORMAT:BINARY";
	else *out_ << "\tFORMAT:TEXT";
	if (concatenated_) *out_ << "\tCONCATENATED";
	*out_ << "\tINDEXED\n";
	*out_ << data.header();
	table_open_=true;
}

template <class T>
void Indexed_file<T>::set_index(const File_index &index)
{
	file_index_=index;
}

template <class T>
File_index Indexed_file<T>::get_index(void) const
{
	return file_index_;
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


template class Data_file <Allele>;
template class Data_file <Population>;
template class Data_file <Locus>;
template class Data_file <Pooled_data>;
template class Data_file <Sample_gof>;
template class Data_file <File_index>;
template class Data_file <Relatedness>;
template class Data_file <Sample_name>;

template class Indexed_file <Allele>;
template class Indexed_file <Population>;
template class Indexed_file <Locus>;
template class Indexed_file <Pooled_data>;

template class Flat_file <Linkage>;
template class Data_file <Linkage>;
template class Indexed_file <Linkage>;

template class Flat_file <Allele>;
template class Flat_file <Population>;
template class Flat_file <Locus>;
template class Flat_file <Pooled_data>;
template class Flat_file <Sample_gof>;
template class Flat_file <File_index>;
template class Flat_file <Relatedness>;
template class Flat_file <Sample_name>;

template class Mpileup_file <Locus>;

