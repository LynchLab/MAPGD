#ifndef STREAMTOOLS_H_
#define STREAMTOOLS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

/* \brief tokenizes a sting seperated by delimiter delim
 */
std::vector<std::string> split(std::string s, char delim);
#endif
