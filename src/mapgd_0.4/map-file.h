/* TODO: Implement.
 * 
 */
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
	//! Reterns a pointer to a new instance of the derived data class.
	Data *read_header(void);	

	//! Reads from the istream.
	/*   
	 *
	 */
	Base_file& read(Data *);	
	Base_file& read(File_index &, Data *);	

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
	 * Returns 0 if unkown. Note, the number of rows in the file will in 
	 * general not be equal to the final position in the file_index, since
	 * many position will be skipped. .
	 */
	size_t size(void) const; 
	bool eof(void);

};

template <class T>
class Data_file : public Base_file {
private :
protected :
	void read_binary(T&);		//!< Read file in binary mode.		
	void write_binary(const T&);	//!< Write in binary mode.

	virtual void read_text(T&){};
	virtual void write_text(const T&){};

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
	Data_file<T>(){
		open_=false;
		table_open_=false;
		read_=false;
		write_=false;
		delim_column_='\t';
		binary_=false;
		filename_="";
		concatenated_=false;
		indexed_=false;
	};			//!< default constructor
	using Base_file::open;

	///! The function that opens a Data_file (if file).
	void open(const char *, const std::ios_base::openmode &);

	///! Doesn't check extnetion.
	void open_extention(const char *, const std::ios_base::openmode &);

//	~Data_file(){};

	///! Appends to a file.
	/*!
	 * If there is a table open, it is imediately closed. 
	 * This is because the information that follows may no 
	 * longer be from the previous class, and the column labels need to 
	 * changed accordingly. If the Flat_file was not opened in 
	 * concatenate to disk, then make a new file with the extention 
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
	Data_file& write(const T &);

	//! Reads a row from the file and advances one row. 
	/*! Returns the class.
	 */
	Data_file& read(T &);				
};

template <class T>
class  Flat_file : public Data_file <T> {
private:
	void read_text(T&);		//!< Read file in text mode.		DONE
	void write_text(const T&);	//!< Write in text mode.	DONE

	using Data_file<T>::out_;	//(const std::ios_base::openmode &);
	using Data_file<T>::in_;	//(const std::ios_base::openmode &);
	using Base_file::write_;	//(const std::ios_base::openmode &);
	using Base_file::binary_;//(const std::ios_base::openmode &);
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
protected:
	File_index file_index_;	//!< A file_index which turns (id0, id1)->pos.

	void read_text(T&);	//!< Read file in text mode.		DONE
	void write_text(const T&);	//!< Write in text mode.	DONE

	using Data_file<T>::out_;	//(const std::ios_base::openmode &);
	using Data_file<T>::in_;	//(const std::ios_base::openmode &);
	using Base_file::write_;	//(const std::ios_base::openmode &);
	using Base_file::binary_;//(const std::ios_base::openmode &);
	using Base_file::table_open_;
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

	id1_t get_pos(const T &) const;
	void write_header(const T&);		//!< Writes a file header.
	T read_header(void);			//!< Reads a file header.
};

/* to provide compatibility with Mpileup_files.*/
template <class T>
class Mpileup_file: public Indexed_file <T> {
private:
	using Indexed_file<T>::file_index_;

	void read_text(T&);	//!< Read file in text mode.			DONE

	using Indexed_file<T>::write_text;	//!< Write in text mode.	DONE
	using Indexed_file<T>::out_;	//(const std::ios_base::openmode &);
	using Indexed_file<T>::in_;	//(const std::ios_base::openmode &);

	using Base_file::open_;
/*	using Base_file::table_open_;
	using Base_file::read_;
	using Base_file::write_;
	using Base_file::delim_column_;
	using Base_file::binary_;
	using Base_file::filename_;
	using Base_file::concatenated_;
	using Base_file::indexed_;*/
public:
/*	Mpileup_file(){
		open_=false;
		table_open_=false;
		read_=false;
		write_=false;
		delim_column_='\t';
		binary_=false;
		filename_="";
		concatenated_=false;
		indexed_=false;
	};			//!< default constructor*/
	using Indexed_file<T>::open;
	using Indexed_file<T>::set_index;	//!< Sets the File_index.		
	using Indexed_file<T>::get_index;	//!< Returns the File_index.

	/*! \brief Returns the position in the file.
	 *
	 * Positions are guarnteed to be unique and increasing for each row. 
	 */

	using Indexed_file<T>::get_pos;
	using Indexed_file<T>::write_header;	//!< Writes a file header.
	T read_header(void);			//!< Reads a file header.
};

//#include "map-file.hxx"

#endif
