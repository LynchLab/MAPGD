#ifndef _DATA_H_
#define _DATA_H_

#include <vector>
#include <string>

class Data {
public:
	Data(){};
	Data(std::vector <std::string> &){};

//	static const std::string file_name;	//!< The default file extention name.	
//      static const std::string table_name;	//!< The destination table in the Db.

	virtual const std::string get_file_name() const {
		fprintf(stderr, "%s:%d: error. No file name define. exiting.\n", __FILE__, __LINE__);
		exit(1);
		return "";
	};	//!< These are neccisary for virutal obj.
        virtual const std::string get_table_name() const{
		fprintf(stderr, "%s:%d: error. No table name define. exiting.\n", __FILE__, __LINE__);
		exit(1);
		return "";
	};

        std::string header(void) const {return "";};

        virtual const std::string sql_header(void) const {return "";};
        virtual const std::string sql_column_names(void) const {return "";};
        virtual const std::string sql_values(void) const {return "";};

        size_t size(void) const {return 0;};        	//!< The size of the ? in bytes.
};

#endif
