#include "map-file.h"

Base_file::Base_file(void){
	open_=false;
	read_=false;
	write_=false;
	delim_column_='\t';
	binary_=false;
	filename_="";
};	

template <class T>
void Data_file<T>::open(const char* filename, const std::ios_base::openmode &mode){
	std::vector<std::string> line;
	line=split_last(filename, '.');
	if (line.back()==T::file_name) filename_=line[0];
	else filename_=filename;
	if (open_){
		std::cerr << __FILE__ << ":" << __LINE__ << ": " << typeid(this).name() << " is already open." << std::endl;
		exit(0);
	}
	if ( mode & std::fstream::in ){
		in_file_.open( (filename_+T::file_name).c_str(), std::ifstream::in);
		if (!in_file_.is_open() ){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open " << filename << T::file_name << " for reading (1)." << std::endl;				
			exit(0);
		};
		in_=&in_file_;
		read_=true;
		openmode_=mode;
	} else if ( mode & std::fstream::out ){
		out_file_.open( (filename_+T::file_name).c_str(), std::ofstream::out);
		if (!out_file_.is_open() ){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open " << filename << T::file_name << " for writing." << std::endl;				
			exit(0);
		};
		out_=&out_file_;
		write_=true;
		openmode_=mode;
	};
	if (mode & std::fstream::binary) binary_=true;
	open_=true;
}

void Base_file::open_no_extention(const char* filename, const std::ios_base::openmode &mode){
	filename_=filename;
	if (open_){
		std::cerr << __FILE__ << ":" << __LINE__ << ": " << typeid(this).name() << " is already open." << std::endl;
		exit(0);
	}
	if ( mode & std::fstream::in ){
		in_file_.open( (filename_).c_str(), std::ifstream::in);
		if (!in_file_.is_open() ){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open " << filename << " for reading (1)." << std::endl;				
			exit(0);
		};
		in_=&in_file_;
		read_=true;
		openmode_=mode;
	} else if ( mode & std::fstream::out ){
		out_file_.open( (filename_).c_str(), std::ofstream::out);
		if (!out_file_.is_open() ){
			std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open " << filename << " for writing." << std::endl;				
			exit(0);
		};
		out_=&out_file_;
		write_=true;
		openmode_=mode;
	};
	if (mode & std::fstream::binary) binary_=true;
	filename_="";
	open_=true;
}

void Base_file::open(const std::ios_base::openmode &mode)
{
	if ( open_ ){
		std::cerr << __FILE__ << ":" << __LINE__ << ": " << typeid(this).name() << " is already open." << std::endl;
		exit(0);
	}
	if ( mode & std::fstream::in ){
		in_=&std::cin;
		read_=true;
		openmode_=mode;
	} else if ( mode & std::fstream::out) {
		out_=&std::cout;
		write_=true;
		openmode_=mode;
	}
	if ( mode & std::fstream::binary) binary_=true;
	open_=true;
}

template <class T>
void Data_file<T>::open_header(Base_file &file){
	if (file.openmode() & std::fstream::in){
		if(file.filename().size()!=0) this->open(file.filename().c_str(), std::fstream::in);
		else this->open(std::fstream::in);
	} else if (file.openmode() & std::fstream::out) {
		if(file.filename().size()!=0) this->open(file.filename().c_str(), std::fstream::out);
		else this->open(std::fstream::out);
	}
	binary_=(file.openmode() & std::fstream::binary);
	open_=true;
}

template <class T>
void Data_file<T>::open_append(Base_file &file)
{
	if(file.is_open() ) file.close();
	if (file.openmode() & std::fstream::in){
		if(file.filename().size()!=0) this->open(file.filename().c_str(), std::fstream::in);
		else this->open(std::fstream::in);
	} else if (file.openmode() & std::fstream::out) {
		if(file.filename().size()!=0) this->open(file.filename().c_str(), std::fstream::out);
		else this->open(std::fstream::out);
	}
	binary_=(file.openmode() & std::fstream::binary);
	open_=true;
}

const std::fstream::openmode& Base_file::openmode(void)
{
	return openmode_;
}

const std::string& Base_file::filename(void)
{
	return filename_;
}

void Base_file::close_table(void){
	if (write_ && open_) *out_ << "@END_TABLE\n";
}

void Base_file::close(void)
{
	if (open_){
		close_table();
		if (write_) {
			out_->setstate(std::ios::eofbit);
			write_=false;
		}
		if (read_){
			in_->setstate(std::ios::eofbit);
			read_=false;
		}
	}
	open_=false;
}

template <class T>
id1_t Indexed_file<T>::get_pos(const T &data) const 
{
	return file_index_.get_rowid(data.id0, data.id1);
}

template <class T>
std::istream& Data_file<T>::read(T &data)
{
	if (!open_) return *in_;
	if (read_){
		if (binary_) read_binary(data);
		else read_text(data);
	} else {
		std::cerr << __FILE__<< ":" <<__LINE__ << ": file not open for reading. The methods Flat_file<type>::open() and Flat_file<type>::read_header(<type>) should be called.";
		exit(0);
	}
//	if (!in_->good() ) { std::cerr << __FILE__ << ":" << __LINE__ << ": unexpected error reading file. Exiting.\n"; exit(0);};
	return *in_;
}

template <class T>
void Flat_file<T>::read_text(T &data)
{
	if (in_->peek()=='@') {
		std::string line;
		*in_ >> line;
		if (line!="@END_TABLE") {
			std::cerr << __FILE__ << ":" << __LINE__ << ": file not closed correctly, exiting.\n";
			exit(0);
		}
		this->close();
	}
	else *in_ >> data;
}

template <class T>
void Indexed_file<T>::read_text(T &data)
{
	std::string scaffold;
	if (in_->peek()=='@') {
		std::string line;
		*in_ >> line;
		if (line!="@END_TABLE") {
			std::cerr << __FILE__ << ":" << __LINE__ << ": file not closed correctly, exiting.\n";
			exit(0);
		}
		this->close();
	} else {
		*in_ >> scaffold;
		if (scaffold.size()==0) in_->setstate(std::ios::failbit);
		if (!(in_->good() ) ) {
			this->close();
		} else {
			data.id0=file_index_.get_id0(scaffold);
			*in_ >> data;
		}
	}
}

template <class T>
void Mpileup_file<T>::read_text(T &data)
{
	std::string scaffold;
	*in_ >> scaffold;
	if (!in_->good() ) {
		this->close();
	} else {
		data.id0=file_index_.get_id0(scaffold);
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
	*out_ << file_index_.get_string(data.id0) << '\t' << data << std::endl;
}


template <class T>
void Data_file<T>::write_binary(const T &data)
{
	out_->write( (const char *) &data, data.size() );
}

template <class T>
std::ostream& Data_file<T>::write(const T &data)
{
	if (binary_) write_binary(data);
	else write_text(data);
//	if (!out_->good() ) { std::cerr << __FILE__ << ":" << __LINE__ << ": unexpected error writing file. Exiting.\n"; exit(0);};
	return *out_;
}

bool Base_file::is_open(void) const{
	return open_;
}			

template <class T>
void Flat_file<T>::write_header(const T &data){
	if (write_) {
		*out_ << "@NAME:" << T::table_name << '\t' <<"VERSION:" << VERSION <<std::endl;
		*out_ << data.header();
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ << ": file not open for writing. Exiting.\n"; exit(0);
	}
}

template <class T>
T Flat_file<T>::read_header(void){
	std::string line;
	std::vector <std::string> columns;
	std::getline(*in_, line);
	columns=split(line, '\t');
	if (columns.size()>0){
		if (columns[0]=="@NAME:"+T::table_name ){
			std::getline(*in_, line);
			columns=split(line, '\t');
			T data(columns);
			return data;
		}
		std::cerr << __FILE__ << ":" << __LINE__ << " attempted to open incorrect header.\n"; 
	}
	std::cerr << __FILE__ << ":" << __LINE__ << " could not initilize " << typeid(T).name() <<"\n";
	T data(columns);
	return data;
}

template <class T>
T Indexed_file<T>::read_header(void)
{
	Flat_file <File_index> index;
	index.open_header(*this);
	file_index_=index.read_header();
	index.read(file_index_);
	index.close_table();
	std::string line;
	std::vector <std::string> columns;
	std::getline(*in_, line);
	columns=split(line, '\t');
	if (columns.size()>0){
		if (columns[0]=="@NAME:"+T::table_name){
			std::getline(*in_, line);
			columns=split(line, '\t');
			T data(columns);
			return data;
		}
		std::cerr << __FILE__ << ":" << __LINE__ << " attempted to open incorrect header.\n"; 
	}
	std::cerr << __FILE__ << ":" << __LINE__ << " could not initilize " << typeid(T).name() << "\n";
	T data(columns);
	return data;
}

template <class T>
T Mpileup_file<T>::read_header(void)
{
	std::string line;
	std::vector <std::string> columns;
	std::getline(*in_, line);
	columns=split(line, '\t');
	T data( (columns.size()-3.)/3. );
//	T data( (columns.size()-2.)/3. );
	return data;
}

template <class T>
void Indexed_file<T>::write_header(const T &data){
	if (!write_) {
		std::cerr << __FILE__ << ":" << __LINE__ << " file not open for writing. Exiting \n";
		exit(0);
	}
	Flat_file <File_index> index;
	index.open_header(*this);
	index.write_header(file_index_);
	index.write(file_index_);
	index.close_table();
	*out_ << "@NAME:" << T::table_name << '\t' <<"VERSION:" << VERSION <<std::endl;
	*out_ << data.header();
}

template <class T>
void Indexed_file<T>::set_index(const File_index &index){
	file_index_=index;
}

template <class T>
File_index Indexed_file<T>::get_index(void) const{
	return file_index_;
}

bool Base_file::eof(void){
	if (open_){
		if (read_) return in_->eof();
		else if (write_) return out_->eof();
		return true;
	}
	else return true;
}

template class Data_file <allele_stat>;
template class Data_file <population_genotypes>;
template class Data_file <Locus>;
template class Data_file <Clone_gof>;
template class Data_file <File_index>;

template class Indexed_file <allele_stat>;
template class Indexed_file <population_genotypes>;
template class Indexed_file <Locus>;

template class Flat_file <allele_stat>;
template class Flat_file <population_genotypes>;
template class Flat_file <Locus>;
template class Flat_file <Clone_gof>;
template class Flat_file <File_index>;

template class Mpileup_file <Locus>;
