#ifndef _DATA_H_
#define _DATA_H_

#include <vector>
#include <string>

class Data {
public:
	Data(std::vector <std::string> &);
	static const std::string file_name;	
        static const std::string table_name;	//!< The destination table in the Db.

        std::string header(void) const {return "";};
        std::string sql_header(void) const {return "";};
        std::string sql_values(void) const {return "";};

        size_t size(void) const {return 0;};        	//!< The size of the ? in bytes.
};

#endif
