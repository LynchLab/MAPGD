#ifndef _DATA_H_
#define _DATA_H_

#include <vector>
#include <string>

class Data {
public:
	Data(){};
	Data(std::vector <std::string> &){};
        virtual const std::string header(void) const {return "";};
	virtual const std::string get_file_name() const {return "";};
        virtual const std::string get_table_name() const{return "";};
        virtual const std::string sql_header(void) const {return "";};
        virtual const std::string sql_column_names(void) const {return "";};
        virtual const std::string sql_values(void) const {return "";};

        size_t size(void) const {return 0;};        	//!< The size of the ? in bytes.
};

#endif
