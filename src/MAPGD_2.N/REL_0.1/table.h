#ifndef _TABLE_H_
#define _TABLE_H_

/*reads and writes tables. These tables should be printible in a csv format (for R) or */

#include <fstream>
#include <cstdio>
#include <iostream>
#include <list>
#include "genotype.h"
#include <map>        

typedef hash_t unsigned int;

class RowIndex {
private:
	std::map <hash_t, TableRow *> _index;
public:
	hash_t (*get_inner)(void);
	hash_t (*get_outer)(void);
	TableRow *get(const RowIndex &);
}

template <class type> class TabelColumn{
private:
	std::string _name;
	size_t _size;
	type *_value;
public:
	int TableColumn(std::string);
	int set(value);
}

class TableHeader {
private:
	std::vector <TableColumn> _columns;
	size_t _size;
public:
	int add (const TableColumn &);
	int del (std::vector <TableColumn>::iterator);
	TableColumn get(const size_t &);
};

class TableRow {
private:
	void * _value;
public:
	int TableRow (const TableHeader &);
	int read (void);
	int write (void);
};

/*Breif: A class to store indexed data.*/
class Table {
private:
	TableRow *row;			//?
	char delim;			//A delimiter used for printing plaintext.
public:
	delim(const char &); 		//sets the delimiter.
	int open(char *);		//opens the file in mode ?
	int open(char *, char *);	//opens the file named ? in mode ?
	int close(void);		//closes the file.
	TableRow read(void);		//reads a row.
	int write(void);		//writes a row.
	int writeheader(void);		//writes the header.
	int readheader(TableHeader);	//reads the header.
	int addcolumn(char);
	int addcolumn(count_t);
	int addcolumn(float_t);
	int seek(RowIndex);		//seeks to index.
}

#endif 
