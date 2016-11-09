#include "sample_name.h"

const std::string Log::file_name=".log";
const std::string Log::table_name="LOGS";
const bool Log::binary=false;

const Registration Log::registered=Registration(Log::table_name, Log::create);

const std::string 
current_time(const time_t &time) 
{
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&time);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}

Log::Log (const std::string &name, const std::string &message)
{
	name_=name;
	message_=message;
	time_=time(0);
	delim='\t';
}

void 
Log::read (std::istream& in) 
{
	std::string line;
	std::getline(in, line);
	std::vector<std::string> columns=split(line, delim_);
	if (columns.size()==3){
		struct tm tm;
		strptime(column[0], "%Y-%m-%d.%X", &tm);
		time_=mktime(&tm);
		name_=column[1];
		message_=column[2];
	} else {
		std::cerr << __FILE__ << ":" <<__LINE__ << ". Table format error."
		std::cerr << " Expected 3 columns, saw " << columns.size() << std::endl;
		exit(0);
	}
}

void 
Log::write(std::ostream& out) const
{
	out << format_time(time) << delim_ << name_ << delim << message_; 
	out << std::endl; 
}

const std::string 
Log::header(void) const 
{
	return "@TIME\tCOMMAND\tMESSAGE\n";
}

const std::string 
Log::sql_header(void) const {
	return "(TIME varchar(255), COMMAND varchar(255), MESSAGE varchar(255) )";
}

const std::string 
Log::sql_column_names(void) const {
	return "(TIME, COMMAND, MESSAGE)";
}

const std::string 
Log::sql_values(void) const {
	return "";
}

size_t 
Log::size(void) const 
{
	return 0;
}

const std::string 
Log::get_file_name(void) const
{
	return file_name;
}

const std::string 
Log::get_table_name(void) const
{
	return table_name;
}

const bool 
Log::get_binary(void) const
{
	return binary;
}
