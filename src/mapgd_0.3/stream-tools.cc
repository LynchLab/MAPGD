#include "stream-tools.h"

std::vector<std::string> split(std::string s, char delim) {
	std::stringstream ss;
	ss << s;
	std::vector<std::string> elems;
	std::string item;
	while (std::getline(ss, item, delim) ) {
		elems.push_back(item);
	//	item.clear();
	}
	return elems;
}
