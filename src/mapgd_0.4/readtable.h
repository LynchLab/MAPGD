#ifndef _READTABLE_HPP_
#define _READTABLE_HPP_

#include <list>
#include "genotype.h"
#include "stream-tools.h"
#include <map>        

typedef std::vector<POPGL> table_t;

table_t readtable(const char *filename);
table_t readbin(const char *filename);
std::map <PAIRGL, size_t> readcounts(const char *filename, size_t a, size_t b, size_t s);
void writebin(table_t table, const char *filename);
void streamtable(const char *, const char *);
#endif 
