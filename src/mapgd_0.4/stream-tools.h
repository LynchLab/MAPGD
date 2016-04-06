#ifndef _STREAMTOOLS_H_
#define _STREAMTOOLS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdio.h>

#ifdef POSIX
#include <sys/poll.h>
#endif 

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

/// An overloaded decleration of split.
inline std::vector<std::string> split(const std::string &s, const char &delim)
{
	std::stringstream ss;
	ss << s;
	return split(ss, delim);
}

///Default behavior is to split on white space and remove it
/*inline std::vector<std::string> split(std::istream &in)
{
	std::string line_in;
	std::string line_out;
	getline(in, line_in);
	line_out.reserve(line_in.size() );
	//TODO REMOVE WHITESPACE
	return split(line_out, ' ');
}*/

/// Takes a string and returns a vector with two elements, split by the delimiter.
/* For example, if the string "A,B,C , D" is given, and a delimiter ',' is used
 * to split the string, then {"A", "B,C , D"} is returned.
 */
inline std::vector<std::string> split_first(std::istream &in, const char &delim)
{
	std::string first, second;
	std::getline(in, first, delim);
	std::getline(in, second);
	if (second.size()!=0) return std::vector <std::string> {first, second};
	else return std::vector <std::string> {first};
}


/// An overloaded decleartion of split_first.
inline std::vector<std::string> split_first(const std::string &s, const char &delim)
{
	std::stringstream ss;
	ss << s;
	return split_first(ss, delim);
}

/// Takes a string and returns a vector with two elements, split by the delimiter.
/* For example, if the string "A,B,C , D" is given, and a delimiter ',' is used
 * to split the string, then {"A,B,C "," D"} is returned.
 */
inline std::vector<std::string> split_last(std::istream &in, const char &delim)
{
	std::string line;
	getline(in, line);
	size_t pos=line.rfind(delim);
	size_t end=line.size();
	if (pos!= std::string::npos ) return std::vector <std::string> {line.substr(0, pos), line.substr(pos, end)};
	else return std::vector <std::string> {line};

}

/// An overloaded decleration of split_last.
inline std::vector<std::string> split_last(const std::string &s, const char &delim)
{
	std::stringstream ss;
	ss << s;
	return split_last(ss, delim);
}

///Removes special characters from a user generated string.
inline std::string sanitize (std::string &s){
	std::replace( s.begin(), s.end(), '\'', '"' );
	std::replace( s.begin(), s.end(), '>', ' ' );
	std::replace( s.begin(), s.end(), '<', ' ' );
	return s;
}

///Removes special characters from a user generated string.
inline std::string sanitize (const std::string &s){
	std::string r=s;
	std::replace( r.begin(), r.end(), '\'', '"' );
	std::replace( r.begin(), r.end(), '>', ' ' );
	std::replace( r.begin(), r.end(), '<', ' ' );
	return r;
}

#ifdef POSIX
///Test a stream to see if it is open
bool check_stream(std::istream *);
#endif

#endif
