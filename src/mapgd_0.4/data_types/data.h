#ifndef _DATA_H_
#define _DATA_H_

#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "typedef.h"

//! A class which can be written as flat text file or into an SQL database.
/*! Data can be written in a plain text representation (the overloaded 
 * \>\> and \<\<), in a binary representation which requires accurate 
 * information from the size() function, or be given to an SQL database 
 * (the sql_ functions). Additionally, all Data must have a static 
 * Registration, which uses the static table_name to write the derived class 
 * into new_data_. The static member new_data_ is a map between strings and 
 * constructors for the derived classes, so that derived Data classes 
 * can be created from strings. Ultimately our goal is to define a Bijection 
 * between the sets {X, Y, ... } and {Z} meaning that we can construct
 * X, Y, ... whenever we have Z and visa versa. However, we cannot know
 * which Data will be necessary to construct Z until X, Y, and Z are all 
 * defined. So, we will define functions whose sets ... injective and surjective function ... Additionally we cannot ...
 * so if {X, Y} is bijective with {Z}, then {X} and {Y} must be surjective on {Z},
 * and {Z} must be injective on both {X} and {Y}. Additionally, if {X, Y} is 
 * surjective on {Z}, we can define a bijective by requiring that ....
 * Additionally, some Data can be transformed into
 * other Data. If the set of {X, Y, ...} is surjective on Z, meaning that 
 * element of Z can be determined from the set {X, Y, ...} then we can 
 * construct Z from the set {X, Y, ...}. If X is injective on Z, then there 
 * exist a unique element of Z for every element of X, such that we  is 
 */
class Data {
private:
	static const std::string file_name;
	static const std::string table_name;
protected:
	//! The read function must be defined in the child class.
	virtual void read(std::istream& str) = 0;
	//! The write function must be defined in the child class.
	virtual void write(std::ostream& str) const = 0;

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
	Data(){};
	Data(std::vector <std::string> &){};
	virtual const std::string get_file_name() const {return this->file_name;};
        virtual const std::string get_table_name() const {return this->table_name;};

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
	//! Return the values to be placed in columns
	/* E.g. 1, 0, ...
 	 */ 
	virtual void sql_read(std::istream &);

	//! Use the \<\< operator to write Data in text mode.
	friend std::ostream& operator << (std::ostream&, const Data &);	
	//! Use the \>\> operator to read Data in text mode.
	friend std::istream& operator >> (std::istream&, Data &);		
	//! The size of the class in bytes. 
	/* If size() is undefined in the child class, then we don't want to 
	 * read or write any at all, so we return 0.
	 */
        virtual size_t size(void) const {return 0;};
	//! Constructs an instance of the class Registered w/ string.
	static Data * new_from_str (const std::string &, const std::vector<std::string> &);
};

//! Data which has an absolute position.
/* Indexed_data can be stored in indexed files.
 */
class Indexed_data : public virtual Data {
protected:
	id1_t abs_pos_;
public:
	Indexed_data(){};
	Indexed_data(std::vector <std::string> &){};
	id1_t get_abs_pos (void) const;
	void set_abs_pos (const id1_t &);
};

//! A static initializer for every translation unit.
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
