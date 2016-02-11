#include "data.h"

//! A 
static int nifty_counter=0; 

std::map <std::string, Data*(*)(const std::vector <std::string> &) > m_data_ctor;
/*static typename std::aligned_storage<sizeof (std::map <std::string, Data*(*)(const std::vector <std::string> Stream) >), alignof (std::map <std::string, Data*(*)(const std::vector <std::string> Stream) >)>::type
  m_data_ctor_buf; // memory for the stream object
std::map <std::string, Data*(*)(const std::vector <std::string> Stream) >& stream = reinterpret_cast<std::map <std::string, Data*(*)(const std::vector <std::string>) >& > (m_data_ctor_buf);*/

Registry_initalizer::Registry_initalizer()
{
	if (nifty_counter++ == 0) new (&m_data_ctor) std::map <std::string, Data*(*)(const std::vector <std::string> &)>; 
}

Registry_initalizer::~Registry_initalizer()
{
	if(nifty_counter-- == 0) free(&m_data_ctor);
}

Data* Data::new_from_str(const std::string &Name, const std::vector<std::string> &columns)
{
	return m_data_ctor[Name](columns);
}

std::istream& operator >> (std::istream& in, Data &data)
{
	data.read(in);
	return in;
}

std::ostream& operator << (std::ostream& out, const Data &data)
{
	data.write(out);
	return out;
}

Registration::Registration (const std::string &str, Data*(*fn)(const std::vector <std::string> &) )
{
#ifdef DEBUG
	std::cerr << "Registering " << str << std::endl;
#endif
	m_data_ctor[str]=fn;
	name_=str;
}


Registration::~Registration ()
{
#ifdef DEBUG
	std::cerr << "Deregistering " << name_ << std::endl;
#endif
	m_data_ctor.erase(name_);
}


id1_t Indexed_data::get_abs_pos (void) const
{
	return abs_pos_;
}

void Indexed_data::set_abs_pos (const id1_t &new_pos)
{
	abs_pos_=new_pos;
}
