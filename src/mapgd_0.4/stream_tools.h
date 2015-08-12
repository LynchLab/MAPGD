#ifndef STREAMTOOLS_H_
#define STREAMTOOLS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdarg>
#include <sstream>

/// Takes a stream/string and returns a vector with one element for each field.
/* 
 * For example, if the string "A,B,C , D" is given, then {"A", "B", "C "," D"} is returned.
 * If the string "A,B,C,\nD,E,F" is given, then {"A", "B", "C"} is returned.
 */
inline std::vector<std::string> split(std::istream &in, const char &delim)
{
	std::vector<std::string> elements;
	std::string line;
	size_t last_delim=0, this_delim=0;
	getline(in, line);
	do { 
		this_delim=line.find(delim, last_delim);
		elements.push_back( line.substr(last_delim, this_delim-last_delim) );
		last_delim=this_delim+1;
	} while ( this_delim!=std::string::npos );
	return elements;
}

inline std::vector<std::string> split(const std::string &s, const char &delim)
{
	std::stringstream ss;
	ss << s;
	return split(ss, delim);
}

inline std::vector<std::string> split(std::istream &in)
{
	std::vector<std::string> elements;
	std::string line;
	size_t last_delim=0, this_delim=0;
	getline(in, line);
	do {
		this_delim=line.find(" \t", last_delim);
		if (this_delim-last_delim>0) elements.push_back( line.substr(last_delim, this_delim-last_delim) );
		last_delim=this_delim+1;
    	} while ( this_delim!=std::string::npos );
	return elements;
}

inline std::vector<std::string> split(const std::string &s)
{
	std::stringstream ss;
	ss <<s;
	return split(ss);
}
	
///Default behavior is to split on white space and remove it

/// Takes a string and returns a vector with two elements, split on first. 
/* For example, if the string "A,B,C , D" is given, then {"A", "B,C , D"} is returned.
 */

/// An overloaded decleartion of split_first.
inline std::vector<std::string> split_first(std::istream &in, const char &delim)
{
	std::string first, second;
	std::getline(in, first, delim);
	std::getline(in, second); 
	return std::vector <std::string> {first, second};
}

inline std::vector<std::string> split_first(const std::string &s, const char &delim)
{
	std::stringstream ss;
	ss << s;
	return split_first(ss, delim);
}

///Stolen from Erik Aronesty's stackoverflow answer.
std::string strprintf(const std::string fmt, ...) 
{
	int size = ((int)fmt.size()) * 2 + 50;   // Use a rubric appropriate for your code
	std::string str;
	va_list args;
	while (true) {     // Maximum two passes on a POSIX system...
		str.resize(size);
		va_start(args, fmt);
		int n = vsnprintf((char *)str.data(), size, fmt.c_str(), args);
		va_end(args);
		if (n > -1 && n < size) {  // Everything worked
			str.resize(n);
			return str;
		}	
		if (n > -1)  // Needed size returned
			size = n + 1;   // For null char
		else
			size *= 2;      // Guess at a larger size (OS specific)
	}
	return str;
}

/// Writes a vector out in a format . . . 
/* For example . . .
 */ 
/*
inline int write(iterator <const std:string> &key, <const std::string> &str, const char &delim)
{
	//SANITIZE STR!!
	out << key << delim << str;
}*/
#endif
