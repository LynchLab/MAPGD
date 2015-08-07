#include "key.h"

KEYNAME::KEYNAME(void)
{
	KEYNAME("NONE");
}

KEYNAME::KEYNAME(const char *this_keyname)
{
	strncpy(keyname_, this_keyname, 8);
	key_=new key("KEYNAME", 0, "a unique lable for the data stored in columns", this);
}

KEYNUM::KEYNUM(void)
{
	KEYNUM(-1);
}

KEYNUM::KEYNUM(uint8_t keynum){
	keyname_="KEYNUM";
	keynum_=1
	keydesc_="a number that is synonymous for a key";
}


KEYDESC::KEYDESC(void)
{
	KEYDESC("NONE");
}

KEYDESC::KEYDESC(std::string desc)
{
	keydesc_=desc;
	key_=new key("KEYDESC", 2, "a verbal description of the data stored by the key", this);
}

///The constructor to be used when keys are declared within data classes.
key::key(const std::string &name, const uint8_t &num, const std::string &desc, data *this_data)
{
	keyname_=KEYNAME(name.c_str());
	keynum_=KEYNUM(num);
	keydesc_=KEYDESC(desc);
	type=&typeid(*this_data);
	size_=this_data->sizeb();
	data_=this_data;
	instance=0;
	offset_=nokey;
}

const char * data::get_name(void) const
{
	return keyname_;
}
const uint8_t data::get_num(void) const
{
	return keynum_;
}
const std::string data::get_desc(void) const
{
	return keydesc_;
}
