#ifndef _DATA_H_
#define _DATA_H_

#include <vector>
#include <string>
#include <map>
#include <iostream>

//! A class which can be written as flat text file or into an SQL database.
/*! Data can be written in a plain text representation (the overloaded 
 * >> and <<) in a binary representation (which requires accurate 
 * information from the size() function, or be given to SQL (the sql_ 
 * functions). Additionally, all Data must have a static Registration, 
 * which uses the static table_name to write the derived class into 
 * new_data_. The static member new_data_ is a map between strings and 
 * constructors for the derived classes, so that derived Data classes 
 * can be created from strings. 
 *
 */
class Data {
protected:
	///! the read function must be ? by the child class.
	virtual void read(std::istream& str)  = 0;
	///! the write function must be ? by the child class.
	virtual void write(std::ostream& str) const = 0;
public:
	Data(){};
	Data(std::vector <std::string> &){};
//      virtual const std::string header(void) const {return "";};
	virtual const std::string get_file_name() const {return "";};
        virtual const std::string get_table_name() const{return "";};
        virtual const std::string sql_header(void) const {return "";};
        virtual const std::string sql_column_names(void) const {return "";};
	virtual const std::string sql_values(void) const {return "";};
	///! use the << operator to write Data in text mode.
	friend std::ostream& operator << (std::ostream&, const Data &);	
	///! use the >> operator to read Data in text mode.
	friend std::istream& operator >> (std::istream&, Data &);		
	//!< The size of the class in bytes.
        virtual size_t size(void) const {return 0;};        	
	static Data * new_from_str (const std::string &, const std::vector<std::string> &);
};

//! A static initilizer for every translation unit.
/*  This initializes the m_data_ctor map which should not be visible
 *  to 
 */
static struct Registry_initalizer {
	Registry_initalizer ();
	~Registry_initalizer ();
} Registry_initalizer;

//! A class which registers a child of Data in Data::new_data_.
/*! All children of Data should have a static declaration of  
 *  Registration so that they can be created by sql_read and 
 *  sql_write.
 */
class Registration 
{
public:
	//!< The constructor of Registration, which stores a constructor for a Data class to new_data_
	Registration (const std::string &str, Data*(*fn)(const std::vector <std::string> &) );
	//!< The constructor of Registration, which stores a constructor for a Data class to new_data_
	~Registration (void);
private:
	std::string name_;
};
#endif
