#ifndef _DATA_H_
#define _DATA_H_

#include <vector>
#include <string>
#include <map>
#include <libintl.h>
#include <iostream>

#include "typedef.h"
#include "error_codes.h"

//! A class which can be written as flat text file or into an SQL database.
/*! Data can be written in a plain text representation (the overloaded 
 * \>\> and \<\<), in a binary representation which requires accurate 
 * information from the size() function, or be given to an SQL database 
 * (the sql_ functions). Additionally, all Data must have a static 
 * Registration, which uses the static table_name to write the derived class 
 * into new_data_. The static member new_data_ is a map between strings and 
 * constructors for the derived classes, so that derived Data classes 
 * can be created from strings. Ultimately my goal is to require Data classes
 * to define a transformations between each other, so that mapgd can convert 
 * between exteranal data schemes and internal schemes as easily as possible, 
 * and maybe even aid in determining when data in the SQL data base needs to 
 * be updated. Notes to self: I can implement his be defining a bijection 
 * between the sets {X, Y, ... } and {Z}, meaning that I can construct 
 * {X, Y, ... } whenever I know {Z} and visa versa. In order to know how to 
 * define a bijection we need to define injective and surjective function 
 * between Data. If {X, Y} are surjective on {Z} but not injective, then we can
 * define a bijection by removing elements from {X} or {Y}. TODO: Find a way to 
 * establish which elements of {X} or {Y} are extraneous.
 * Also, is this injective/surjective language actually helpful?
 */
class Data {
private:
	/*! A flag to indicate whether reading/writing in binary mode is 
	 * supported for the data_type.
 	 */
        static const bool binary;
protected:

	//! The read function must be defined in the child class.
	virtual void read(std::istream& str) = 0;

	//! The write function must be defined in the child class.
	 virtual void write(std::ostream& str) const = 0;

//	Data
/*	//! T can be constructed form this, but this can not be constructed from T.
	template <class T, typename... Args> 
	void 
	injection_ (T&) const { 
		injection_{}
	}
	//! T cannot be constructed form this, but this can be constructed from T.
	template <class T> 
	void
	surjection_ (const T&) { }

	template<typename T, typename... Args> 
	T& convert(const Args &args ...)
	{
		
		args.injection_(T);
	}*/

public:
	//! The default file extension to use.
	static const std::string file_name;
	//! The name of the table in SQL databases.
	static const std::string table_name;

	//! The list of any deprecated table names.
	static const std::vector <std::string> table_names_old;

	virtual std::string header(void) const = 0;

	void read_binary(std::istream& str) {};
	void write_binary(std::ostream& str) const {};

	Data(){};
	Data(std::vector <std::string> &){};

	virtual const std::string get_file_name() const {return this->file_name;};
        virtual const std::string get_table_name() const {return this->table_name;};

	virtual const bool indexed() const {return false;};

        virtual const bool get_binary() const {return this->binary;};

	/**
	 * \defgroup SQL functions
	 * @{
	 */

	//! Return the names of the columns, along with variable type
	/* E.g. Column_name_1 integer, Column_name_2 integer, ...
 	 */ 
        virtual const std::string sql_header(void) const {return "";};
	//! Return the names of the columns
	/* E.g. Column_name_1, Column_name_2, ...
 	 */ 
        virtual const std::string sql_column_names(void) const {return "";};
	//! Return the values to be placed in columns
	/* E.g. 1, 0, ...
 	 */ 
	virtual const std::string sql_values(void) const {return "";};
	//! Reads the values...
	/* E.g. 1, 0, ...
 	 */ 
	virtual void sql_read(std::istream &);

	/*
	 *@}
         */

	/**
	 * \defgroup basic IO
	 * @{
	 */
	//! Use the \<\< operator to write Data in text mode.
	friend std::ostream& operator << (std::ostream&, const Data &);	
	//! Use the \>\> operator to read Data in text mode.
	friend std::istream& operator >> (std::istream&, Data &);		

	//! The size of the class in bytes. 
	/* If size() is undefined in the child class, then we don't want to 
	 * read or write any at all, so we return 0.
	 */
        virtual size_t size(void) const {return 0;};
	/*
	 *@}
         */

	//! Constructs an instance of the class Registered w/ string.
	static Data * new_from_str (const std::string &, const std::vector<std::string> &);
//	static Data * new_empty (const std::string &);

	//! Toys
//	virtual const enum tags(void) const {return "";};

};

//! Data which has an absolute position.
/* Indexed_data can be stored in indexed files.
 */
class Indexed_data : public virtual Data {
protected:
	id1_t abs_pos_;
public:
	using Data::write_binary;
	using Data::read_binary;
	void write_pos(std::ostream& str) const
	{
		str.write((char *)(&abs_pos_), sizeof(id1_t) );
	}

	void read_pos(std::istream& str)
	{
		str.read((char *)(&abs_pos_), sizeof(id1_t) );
	}

	Indexed_data(){};
	Indexed_data(std::vector <std::string> &){};
	id1_t get_abs_pos (void) const;	  //!< Indexed data needs to associate each datum with a position in the genome.

	void set_abs_pos (const id1_t &); //!< Indexed data needs to associate each datum with a position in the genome. 
	const bool indexed() const {return true;};
};


//! Data which has an absolute position.
/* Indexed_data can be stored in indexed files.
 */

#define get_abs_pos1	get_abs_pos
#define set_abs_pos1	set_abs_pos

class Double_indexed_data : public virtual Data 
{
protected:
	id1_t abs_pos1_, abs_pos2_;
public:
	using Data::write_binary;
	using Data::read_binary;

	void write_pos(std::ostream& str) const
	{
		str.write((char *)(&abs_pos1_), sizeof(id1_t) );
		str.write((char *)(&abs_pos2_), sizeof(id1_t) );
	}

	void read_pos(std::istream& str)
	{
		str.read((char *)(&abs_pos1_), sizeof(id1_t) );
		str.read((char *)(&abs_pos2_), sizeof(id1_t) );
	}

	Double_indexed_data(){};
	Double_indexed_data(std::vector <std::string> &){};

	id1_t get_abs_pos (void) const;	  //!< Indexed data needs to associate each datum with a position in the genome.
	id1_t get_abs_pos2 (void) const;	  //!< Indexed data needs to associate each datum with a position in the genome.

	void set_abs_pos (const id1_t &); //!< Indexed data needs to associate each datum with a position in the genome. 
	void set_abs_pos2 (const id1_t &); //!< Indexed data needs to associate each datum with a position in the genome. 

	const bool indexed() const {return true;};
};

//! A static initializer for every translation unit.
/*  This initializes m_data_ctor map which should not be visible
 *  outside of data.cc
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

std::vector <std::string> registry_list(void);

#endif
