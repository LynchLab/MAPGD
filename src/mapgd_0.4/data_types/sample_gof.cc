#include "sample_gof.h"

const std::string Sample_gof::file_name=".gof";
const std::string Sample_gof::table_name="SAMPLE";
const bool Sample_gof::binary=false;

const Registration Sample_gof::registered=Registration(Sample_gof::table_name, Sample_gof::create);

Sample_gof::Sample_gof ()
{
	name_="";
	number_=0;
	smp_num_=0;
	delim='\t';
}

Sample_gof::Sample_gof (const std::string &name, const float_t &number)
{
	name_=name;
	number_=number;
	smp_num_=0;
	delim='\t';
}

Sample_gof::Sample_gof (const count_t &id, const std::string &name, const float_t &number)
{
	//std::cerr << id << ", " << sanitize(name) << std::endl;
	name_=name;
	number_=number;
	smp_num_=id;
	delim='\t';
}

void
Sample_gof::read (std::istream& in)
{
	std::string line, temp;
	std::getline(in, line);
	std::stringstream line_stream(line);
	//std::cerr << "read line:" << line << std::endl;
	line_stream >> smp_num_ >> name_ >> number_;
	if (!line_stream.eof() )
	{
		line_stream.clear();
		std::string temp;
		line_stream >> smp_num_ >> name_ >> temp;

		if(temp!="." )//std::string(MISSING))
		{
			//std::cerr << temp << " - " << bool(temp==".") << std::endl;
			fprintf(stderr, gettext("mapgd:%s:%d: Sample_gof::read: error parsing input\n"), __FILE__, __LINE__ );
			exit(BADIN);
		} else {
			//std::cerr << temp << " - " << bool(temp==".") << std::endl;
		}
	}
	/*
	*/
}

void
Sample_gof::sql_read (std::istream& in)
{
	std::string line, temp;
	std::getline(in, line);
	std::stringstream line_stream(line);
	line_stream >> smp_num_ >> name_ >> number_;
}

void
Sample_gof::write (std::ostream& out) const
{
	if (!isnan(number_) )
		out << smp_num_ << delim << name_ << delim << number_;
	else 
		out << smp_num_ << delim << name_ << delim << MISSING;
}

std::string 
Sample_gof::header(void) const 
{
	return std::string("@SMPNUM\tSMPNAME\tGOF\n");
}

size_t 
Sample_gof::size(void) const 
{
	//Oh, lets just have a segfault cause I'm bored.
	return sizeof(float_t)+sizeof(char)+name_.size()*sizeof(char);
}

const bool Sample_gof::get_binary(void) const
{
	return binary;
}

const std::string Sample_gof::get_file_name(void) const
{
	return file_name;
}

const std::string Sample_gof::get_table_name(void) const
{
	return table_name;
}

const std::string Sample_gof::sql_header(void) const {
	return "(SMPNUM int, SMPNAME VARCHAR(255), GOF REAL, PRIMARY KEY (SMPNUM) )";
}

const std::string Sample_gof::sql_column_names(void) const {
	return "(SMPNUM, SMPNAME, GOF)";
}

const std::string Sample_gof::sql_values(void) const {
        char return_buffer[SQL_LINE_SIZE]={0};
#if FLT_EVAL_METHOD == 2
	snprintf(return_buffer, SQL_LINE_SIZE, "(%d, '%s', %Lf)",
	smp_num_,
	name_.c_str(),
	number_);
#elif FLT_EVAL_METHOD == 1
	snprintf(return_buffer, SQL_LINE_SIZE, "(%d, '%s', %f)",
	smp_num_,
	name_.c_str(),
	number_);
#elif FLT_EVAL_METHOD == 0
	snprintf(return_buffer, SQL_LINE_SIZE, "(%d, '%s', %f)",
	smp_num_,
	name_.c_str(),
	number_);
#endif
        return std::string(return_buffer);
}
