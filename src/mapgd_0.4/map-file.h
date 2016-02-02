/* TODO: Implement.
 * 
 */
#ifndef MAP_FILE_H_
#define MAP_FILE_H_	

#include <map>
#include <string>
#include <iostream>
#include <typeinfo>

#include "stream-tools.h"
#include "typedef.h"
#include "stream-tools.h"
#include "allele_stat.h"
#include "locus.h"
#include "genotype.h"
#include "sample_gof.h"
#include "file-index.h"
#include "relatedness_data.h"
#include "pooled_data.h"

// PLEASE LIMIT LINE LENGTH TO 79 CHARACTERS----------------------------------/

/// A templet which stores data associated with specific locations in a genome.
/* An indexed file is a table with colulmns specified by the data types stored
 * in the table, which must store a pair of IDs retrived by the get_id0() and 
 * get_id1() member functions. Additionally, data types must return a name for 
 * the table in the database where the data may be stored with member function 
 * table_name(), and must name the columns of the table with cannonical data 
 * types listed in the file ?. This allows data to easily be transfered into 
 * and out of a database. Finally...
 *
 */

class Base_file {
protected:
	bool open_;	//!< indicates whether the profile opened succesfully
	bool table_open_;	//!< indicates whether the profile opened succesfully

	char delim_column_;	//!< The delimiter which seperates columns

	bool read_;		//!< File is open for reading.
	bool write_;		//!< File is open for writing.
	bool binary_;		//!< Binary mode flag. 

	std::istream *in_;	//!< All data is read from in.
	std::ostream *out_;	//!< All data is writen is writen to out.

	std::string filename_;
	std::fstream in_file_;	//!< The file to read data from (if not stdin).
	std::ofstream out_file_;//!< The file to write data to (if not stdout).

	std::fstream::openmode openmode_;
public:
	const std::fstream::openmode& openmode();//1(openmode_);
	const std::string& filename();//1(openmode_);

	Base_file();				//!< default constructor

	void open_no_extention(const char *, const std::ios_base::openmode &);//!< The function that opens a indexed_file (if file).
	void open(const std::ios_base::openmode&);//!< The function that opens a indexed_file (if stdin).

	void close(void); //!< Close iostreams, writes tail, etc.
	void close_table(void); //!< Ends table w/o closing iostreams.
	bool is_open(void) const;//!< Returns true iff indexed_file is open.
	/*done*/

	/** @defgroup BasicIO Basic IO operations
	 * @{
  	 */

	std::istream& seekg(id1_t pos);		
	std::ostream& seekp(id1_t pos);

	/*! \brief NOT IMPLEMENTED!!!
	 */
	std::istream& seekg(id1_t off, std::ios_base::seekdir way);	

	/*! \brief NOT IMPLEMENTED!!!
	 */
	std::ostream& seekp(id1_t off, std::ios_base::seekdir way);	

	id1_t tellp(void);		//!< Tells row number of puts. 
	id1_t tellg(void);		//!< Tells row number of gets.

	/**@defgroup Formating Formating options
	 * @{
	 */
	void set_delim (const char&);				//!< Sets the delimiter that seperates columns. Only used in text mode.
	const char & get_delim (const char&) const;		//!< Gest the delimiter that seperates columns. Only used in text mode.
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
	void read_binary(T&);	//!< Read file in binary mode. 		DONE
	void write_binary(const T&);	//!< Write in binary mode.	DONE

	virtual void read_text(T&){};
	virtual void write_text(const T&){};

//	using Base_file::out_file_;//(const std::ios_base::openmode &);
//	using Base_file::in_file_;//(const std::ios_base::openmode &);

	using Base_file::out_;//(const std::ios_base::openmode &);
	using Base_file::in_;//(const std::ios_base::openmode &);

public :
	using Base_file::open;//(const std::ios_base::openmode &);
//	using Base_file::eof;
//	using Base_file::close;
	void open(const char *, const std::ios_base::openmode &);//!< The function that opens a indexed_file (if file).

	/*! \brief Appends to a file.
	 *
	 * If there is a Flat_file open, it is imediately closed. 
	 * This is because the information that follows may no 
	 * longer be from the previous class, and the column labels 
	 * need to changed accordingly.If the Flat_file was opened 
	 * to disk, then make a new file with the extention 
	 * T::file_name is created.
	 */
	void open_append(Base_file &); 

        /*! \brief Opens a header for a Flat_file.
         *
         * This method opens a file to the same i/o stream as 
	 * Base_file iff Base_file is not associated with a file,
	 * and creates a file with the extention T:file_name 
	 * otherwise. Base_file is not closed. This can result in 
	 * data between the Base_file and the . . . being 
	 * interspersed, potentially curupting the file.
         */
        void open_header(Base_file &);

	/*! \brief Writes a row to the file and advances one row. 
	 *
	 * Returns the ostream.
	 */
	std::ostream& write(const T &);
	/*! \brief Reads a row from the file and advances one row. 
	 *
	 * Returns the istream.
	 */
	std::istream& read(T &);				
};

template <class T>
class  Flat_file : public Data_file <T> {
private:
	void read_text(T&);		//!< Read file in text mode.		DONE
	void write_text(const T&);	//!< Write in text mode.	DONE

	using Data_file<T>::out_;	//(const std::ios_base::openmode &);
	using Data_file<T>::in_;	//(const std::ios_base::openmode &);
	using Base_file::write_;	//(const std::ios_base::openmode &);

public:
	using Data_file<T>::open;
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

public:
	using Data_file<T>::open;

	void set_index(const File_index&);	//!< Sets the File_index.		

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

public:
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

#endif
