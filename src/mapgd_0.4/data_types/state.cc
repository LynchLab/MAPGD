#include "state.h"

const std::string State::file_name=".stt";
const std::string State::table_name="STATES";
const bool State::binary=true;

const Registration State::registered=Registration(State::table_name, State::create);
/** @breif constuctor w/ initial values. **/

State::State(const std::vector <std::string> &column_names)
{
	sample_names_=std::vector <std::string> (column_names.cbegin(), column_names.cend() );
	?
}


State::State()
{
	data=NULL;
}

State::State(const State &rhs)
{
	*this=rhs;
}

State::~State()
{
	if (data!=NULL)
		free data;
	data=NULL;
}

/**@breif return size of State if State is set, 0 otherwise**/
size_t 
State::size() const
{
	if (data!=NULL)
	return ?;
	return 0;
}

State & 
State::operator= (const State& rhs)
{
	return *this;
}

std::string 
State::header(void) const
{
	std::vector <std::string>::const_iterator s_it=sample_names_.cbegin(), end=sample_names_.cend();
	while(s_it!=end){
		line+='\t';	
		line+=*s_it;
		s_it++;	
	}
	line+='\n';
	return line;
}

void
State::write (std::ostream& out) const
{
	it_=data;
	while(it_!=end_){
		out << '\t' << get_next();
	}
}

void
State::read (std::istream& in)
{
        std::string line;
        std::getline(in, line);
        std::stringstream line_stream(line);
	trait;
	while(s_it!=end){
		line_stream >> trait;
		set_next(trait);
	}	
}


void
State::write_binary (std::ostream& out) const
{
        out.write((char *)&data, (size_t)(size_*sizeof(char) ) );
}

void
State::read_binary (std::istream& in)
{
        in.read((char *)&data, (size_t)(size_*sizeof(char) ) );
}


const bool
State::get_binary(void) const
{
	return binary;
}


static uint32_t Mask[32]={0x00000001, 0x00000002, 0x00000004, 0x00000008,
                        0x00000010, 0x00000020, 0x00000040, 0x00000080,
                        0x00000100, 0x00000200, 0x00000400, 0x00000800,
                        0x00001000, 0x00002000, 0x00004000, 0x00008000,
                        0x00010000, 0x00020000, 0x00040000, 0x00080000,
                        0x00100000, 0x00200000, 0x00400000, 0x00800000,
                        0x01000000, 0x02000000, 0x04000000, 0x08000000,
                        0x10000000, 0x20000000, 0x40000000, 0x80000000};

uint32_t get(const uint32_t &x) const
{
        if (( (data+x) & Mask[bit_])!=0 ){
                return 1;
        } else  {
                return 0;
        }
}

uint32_t get_next()
{
        if (( (++it_) & Mask[bit_])!=0 ){
                return 1;
        } else  {
                return 0;
        }
}

void set(const uint32_t &x, const uint32_t &v) 
{
        if (( (data+x) & Mask[bit_])!=0 ){
}

void set_next(const uint32_t &v)
{
        if (( (++it_) & Mask[bit_])!=0 ){
}

void set_next32(const uint32_t &x) 
{
        ++it_=x;
}

