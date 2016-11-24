#ifndef _MAP_FILE_H_
#define _MAP_FILE_H_	

#include <map>
#include <string>
#include <iostream>
#include <typeinfo>

#include "stream-tools.h"
#include "typedef.h"
#include "stream-tools.h"
#include "file-index.h"


/* These should all be included from the objects using this code, but I'm not 
 * linking things correctly right now.
 */
#include "allele.h"
#include "linkage_data.h"
#include "locus.h"
#include "genotype.h"
#include "sample_gof.h"
#include "relatedness_data.h"
#include "population.h"
#include "pooled_data.h"
#include "sample_name.h"
#include "bcf2pro.h"

#define READ	std::ios::in
#define WRITE	std::ios::out

#define OPEN	1
// PLEASE LIMIT LINE LENGTH TO 79 CHARACTERS----------------------------------/

//! A templet which stores data associated with specific locations in a genome.
/*! An indexed file is a table with colulmns specified by the data types stored
 * in the table, which must store a pair of IDs retrived by the get_id0() and 
 * get_id1() member functions. Additionally, data types must return a name for 
 * the table in the database where the data may be stored with member function 
 * table_name(), and must name the columns of the table with cannonical data 
 * types listed in the file ?. This allows data to easily be transfered into 
 * and out of a database. Finally...
 *
 */
class Base_file {
private:
protected:

	//! indicates whether the iostream opened succesfully
	bool open_;		

	//! indicates whether the header has been read succesfully
	/*   Since all tables begin with a header, 
	 *   Base_file::read(Data *) will return NULL before until 
	 *   read_header has been called and moved the iostream past 
	 *   the header.
	 */
	bool table_open_;	

	//! indicates whether multiple tables are included in the file
	bool concatenated_;	

	//! The delimiter which seperates columns
	char delim_column_;	

	bool read_;		//!< File is open for reading.
	bool write_;		//!< File is open for writing.
	bool binary_;		//!< Binary mode flag. 
	bool try_binary_;	//!< Attempt to set binary mode flag. 
	bool indexed_;		//!< Indexed mode flag.

	std::istream *in_;	//!< All data is read from in.
	std::ostream *out_;	//!< All data is writen is writen to out.

	std::fstream file_;	//!< The file to read data from.

	//! stores information about the mode in which the file was opened.
	std::ios::openmode openmode_;
	std::string filename_;	//!< The name of the file if opened.
public:
	//! returns the open mode.
	const std::fstream::openmode& openmode();
	//! returns the filename.
	const std::string& filename();		
	//! returns the concatenated flag.
	const bool & concatenated();		

	Base_file();		//!< default constructor

	//! The function that opens a indexed_file (if file).
	void open_no_extention(const char *, const std::ios::openmode &);
	//! Opens a Base_file to the cin/cout.
	void open(const std::ios::openmode&);
	//! Opens a Base_file to an istream.
	void open(std::iostream*, const std::ios_base::openmode &);
	//! Opens a Base_file to an istream.
	void open(std::istream*, const std::ios_base::openmode &);
	//! Opens a Base_file to an istream.
	void open(std::ostream*, const std::ios_base::openmode &);

	//! returns the istream.
	std::istream* get_in(void);
	//! returns the ostream.
	std::ostream* get_out(void);	

	//! Close iostreams, writes tail, etc.
	void close(void); 
	//! Ends table w/o closing iostreams.
	void close_table(void); 
	//! Returns true iff Base_file is open.
	bool is_open(void) const;
	//! Returns true iff a table is open in Base_file.
	bool table_is_open(void) const;
	//! Returns true iff a table is open in binary mode.
	bool binary(void) const;

	/** @defgroup BasicIO Basic IO operations
	 * @{
  	 */
	//! Returns a pointer to a new instance of the derived data class.
	Data *read_header(void);	

	//! Reads from the istream.
	/*   
	 *
	 */
	Base_file& read(Data *);	
	Base_file& read(File_index &, Indexed_data *);	

	//! Writes to the ostream.	
	Base_file& write(Data *);	
	Base_file& write(File_index &, Data *);	

	//! sets in and out (for reading and writing) to possition pos.	
	void seek(const std::streampos &pos);		
	//! sets in (for reading) to possition pos.
	void seekg(const std::streampos &pos);		
	//! sets out (for writing) to possition pos.
	void seekp(const std::streampos &pos);		

	void seek(std::streampos pos, std::ios_base::seekdir way);	//!< TODO		
	void seekg(std::streampos  off, std::ios_base::seekdir way);	//!< TODO
	void seekp(std::streampos  off, std::ios_base::seekdir way);	//!< TODO

	std::streampos tellp(void);	//!< Tells streampos of out (writing) 
	std::streampos tellg(void);	//!< Tells streampos of in (reading)
	/** @} */

	/**@defgroup Formating Formating options
	 * @{
	 */
	//! Sets the delimiter that seperates columns. Only used in text mode.
	void set_delim (const char&);				
	//! Gest the delimiter that seperates columns. Only used in text mode.
	const char & get_delim (const char&) const;		
	/** @} */

	/*functions dealing with ?*/

	/*! \brief Returns the number of rows in the file.
	 *
	 * Returns 0 if unknown. Note, the number of rows in the file will in 
	 * general not be equal to the final position in the file_index, since
	 * many position will be skipped. .
	 */
	size_t size(void) const; 
	bool eof(void);
	bool indexed(void);
	bool check_concatenated(const char *filename);
};

template <class Data>
class Data_file : public Base_file {
private :
protected :
	void read_binary(Data &);		//!< Read file in binary mode.		
	void write_binary(const Data &);	//!< Write in binary mode.

	virtual void read_text(Data &){};
	virtual void write_text(const Data&){};

	using Base_file::out_;
	using Base_file::in_;

	using Base_file::open_;
	using Base_file::table_open_;
	using Base_file::read_;
	using Base_file::write_;
	using Base_file::delim_column_;
	using Base_file::binary_;
	using Base_file::filename_;
	using Base_file::concatenated_;
	using Base_file::indexed_;
public:
	using Base_file::open;

	///! The function that opens a Data_file (if file).
	void open(const char *, const std::ios_base::openmode &);

	///! Doesn't check extension.
	void open_extention(const char *, const std::ios_base::openmode &);

//	~Data_file(){};

	///! Appends to a file.
	/*!
	 * If there is a table open, it is immediately closed. 
	 * This is because the information that follows may no 
	 * longer be from the previous class, and the column labels need to 
	 * changed accordingly. If the Flat_file was not opened in 
	 * concatenate to disk, then make a new file with the extension 
	 * T::file_name is created.
	 */
	void open_from(Base_file &); 

        //! Opens a header for a Flat_file.
        /*! 
         * This method opens a file to the same i/o stream as 
	 * Base_file iff Base_file is not associated with a file,
	 * and creates a file with the extention T:file_name 
	 * otherwise. Base_file is not closed. This can result in 
	 * data between the Base_file and the . . . being 
	 * interspersed, potentially curupting the file.
         */
        void open_header(Base_file &);

	//! Writes a row to the file and advances one row. 
	/*! Returns the ostream.
	 */
	Data_file& write(const Data &);

	//! Reads a row from the file and advances one row. 
	/*! Returns the class.
	 */
	Data_file& read(Data &);				
};

template <class T>
class  Flat_file : public Data_file <T> {
private:
	void read_text(T&);		//!< Read file in text mode.	DONE
	void write_text(const T&);	//!< Write in text mode.	DONE

	using Data_file<T>::out_;	//(const std::ios_base::openmode &);
	using Data_file<T>::in_;	//(const std::ios_base::openmode &);
	using Base_file::write_;	//(const std::ios_base::openmode &);
	using Base_file::binary_;//(const std::ios_base::openmode &);
	using Base_file::try_binary_;//(const std::ios_base::openmode &);
	using Base_file::concatenated_;//(const std::ios_base::openmode &);
	using Base_file::table_open_;

public:
	using Data_file<T>::open;
//	~Flat_file(){};
	void write_header(const T&);		//!< Writes a file header.
	T read_header(void);			//!< Reads a file header.
};
	
template <class T>
class Indexed_file: public Data_file <T> {
private:
	id1_t reference_size_;
protected:
	File_index file_index_;	//!< A file_index which turns (id0, id1)->pos.

	void read_text(T&);	//!< Read file in text mode.		DONE
	void write_text(const T&);	//!< Write in text mode.	DONE

	void read_binary(T &);
	void write_binary(const T &);

	using Data_file<T>::out_;	//(const std::ios_base::openmode &);
	using Data_file<T>::in_;	//(const std::ios_base::openmode &);
	using Base_file::write_;	//(const std::ios_base::openmode &);
	using Base_file::binary_;//(const std::ios_base::openmode &);
	using Base_file::try_binary_;//(const std::ios_base::openmode &);
	using Base_file::table_open_;
	using Base_file::read_;
	using Base_file::concatenated_;//(const std::ios_base::openmode &);
	using Base_file::filename_;

public:
	using Data_file<T>::open;

	void set_index(const File_index&);	//!< Sets the File_index.		
	~Indexed_file(){};

	File_index get_index(void) const;	//!< Returns the File_index.

	/*! \brief Returns the position in the file.
	 *
	 * Positions are guarnteed to be unique and increasing for each row. 
	 */
	Indexed_file& write(const T&);//(const Indexed_data &);
	Indexed_file& read(T&);//(const Indexed_data &);

	id1_t get_pos(const T &) const;
	void write_header(const T&);		//!< Writes a file header.
	T read_header(void);			//!< Reads a file header.
};

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
//		if ( Base_file::check_concatenated(filename) ){
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

template <class T>
void Data_file<T>::open_from(Base_file &file)
{
	if (file.table_is_open() ) file.close_table();
	if (file.openmode() & std::ios::in){
//		if(file.filename().size()!=0) {
		if(file.concatenated() ) this->open(file.get_in(), file.openmode() );
		else this->open_extention(file.filename().c_str(), file.openmode() );
//		}
	} else if (file.openmode() & std::ios::out) {
//		if(file.filename().size()!=0) {
		if(file.concatenated() ) this->open(file.get_out(), file.openmode() );
		else this->open_extention(file.filename().c_str(), file.openmode() );
//		}
	}
	try_binary_=(file.openmode() & std::ios::binary);
	open_=true;
}

template <class T>
id1_t Indexed_file<T>::get_pos(const T &data) const 
{
	return data.get_abs_pos();
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
Indexed_file<T>& Indexed_file<T>::read(T &data)
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
void Data_file<T>::read_binary(T &data)
{
	char a='a';
	if (in_->peek()!='@'){
		in_->read(&a, sizeof(char) );
		data.read_binary(*in_);
	} else {
		if (!concatenated_) {
			this->close();
		} else {
			this->close_table();
		}
	}
}

template <class T>
void Indexed_file<T>::read_binary(T &data)
{
	char a='a';
	if (in_->peek()=='@')
	{
		if (!concatenated_) {
			this->close();
		} else {
			this->close_table();
		}
	} else {
		in_->read(&a, sizeof(char) );
		data.read_pos(*in_);
		data.read_binary(*in_);
	}
	
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
	char a='a';
	out_->write(&a, sizeof(char) );
	data.write_binary(*out_);
}

template <class T>
void Indexed_file<T>::write_binary(const T &data)
{
	char a='a';
	out_->write(&a, sizeof(char) );
	data.write_pos(*out_);
	data.write_binary(*out_);
}

template <class T>
Data_file<T>& Data_file<T>::write(const T &data)
{
	if (binary_) write_binary(data);
	else write_text(data);
//	if (!out_->good() ) { std::cerr << __FILE__ << ":" << __LINE__ << ": unexpected error writing file. Exiting.\n"; exit(0);};
	return *this;
}

template <class T>
Indexed_file<T>& Indexed_file<T>::write(const T &data)
{
	if (binary_) write_binary(data);
	else write_text(data);
//	if (!out_->good() ) { std::cerr << __FILE__ << ":" << __LINE__ << ": unexpected error writing file. Exiting.\n"; exit(0);};
	return *this;
}

template <class T>
void Flat_file<T>::write_header(const T &data)
{
	if (write_) {
		*out_ << "@NAME:" << T::table_name << "\tVERSION:" << VERSION;
		if (try_binary_ && T::binary)
		{
			*out_ << "\tFORMAT:BINARY";
			binary_=true;
		} else {
			*out_ << "\tFORMAT:TEXT";
			binary_=false;
		}
		if (concatenated_) *out_ << "\tCONCATENATED";
		*out_ << std::endl;
		*out_ << data.header();
		table_open_=true;
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ << ": file not open for writing. Exiting.\n"; exit(0);
	}
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
			binary_=std::find(columns.begin(), columns.end(), "FORMAT:BINARY")!=columns.end();
			concatenated_=std::find(columns.begin(), columns.end(), "CONCATENATED")!=columns.end();
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
			binary_=std::find(columns.begin(), columns.end(), "FORMAT:BINARY")!=columns.end();
			concatenated_=std::find(columns.begin(), columns.end(), "CONCATENATED")!=columns.end();
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
	reference_size_=file_index_.get_reference_size();
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
	if (try_binary_ && T::binary){
		*out_ << "\tFORMAT:BINARY";
		binary_=true;
	} else {
		*out_ << "\tFORMAT:TEXT";
		binary_=false;
	}
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

#endif
