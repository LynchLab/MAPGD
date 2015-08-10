#include "key.h"

/* La la La la La la ...  LA la La la La la. SAIL 
 * I'm not really sure this is the right way to do things.
 * I'm thinking it will be nice if keys can be pack into 
 * just like any other data, and makeing these global will
 * probably fundementally cripple my entire implementation, 
 * but, I'm tired.
 */

constexpr uint8_t key::nokey;

KEYNAME::KEYNAME(void)
{
	std::cerr << "this\n";
	KEYNAME("NONE");
}

KEYNAME::KEYNAME(const char *keyname)
{
	std::cerr << "KEYNAME\n";
//	keyname_=&_keyname0;
//	keynum_=&_keynum0;
//	keydesc_=&_keydesc0;
	value_=keyname;
	instance_=0;
	offset_=key::nokey;
	size_=size();
//	size_of_f=&size_of_char7;
//	to_text_f=&to_text_char7;
//	from_tx_f=&from_tx_char7;
}

KEYNUM::KEYNUM(void)
{
	KEYNUM(nokey);
}

KEYNUM::KEYNUM(const uint8_t &keynum)
{
	std::cerr << "KEYNUM\n";
//	keyname_=&_keyname1;
//	keynum_=&_keynum1;
//	keydesc_=&_keydesc1;
	value_=&keynum;
	instance_=0;
	offset_=key::nokey;
	size_=size();
//	size_of_f=&size_of_char7;
//	to_text_f=&to_text_char7;
//	from_tx_f=&from_tx_char7;
}


KEYDESC::KEYDESC(void)
{
	KEYDESC("NONE");
}

KEYDESC::KEYDESC(const char *desc)
{
//        keyname_=&_keyname3;
//        keynum_=&_keynum3;
//        keydesc_=&_keydesc3;
	std::cerr << "KEYDESC\n";
        value_=new std::string(desc);
        instance_=0;
        offset_=nokey;
        size_=0;
//	size_of_f=&size_of_char7;
//	to_text_f=&to_text_char7;
//	from_tx_f=&from_tx_char7;
}

KEYDESC::KEYDESC(const std::string &desc)
{
//        keyname_=&_keyname3;
//        keynum_=&_keynum3;
//        keydesc_=&_keydesc3;
	std::cerr << "KEYDESC\n";
        value_=&desc;
        instance_=0;
        offset_=nokey;
        size_=0;
//	size_of_f=&size_of_char7;
//	to_text_f=&to_text_char7;
//	from_tx_f=&from_tx_char7;
}

key::key()
{
	keyname_=NULL;
	keynum_=NULL;
	keydesc_=NULL;
	value_=NULL;
	instance_=0;
	offset_=key::nokey;
	size_of_f=NULL;
	to_text_f=NULL;
	from_tx_f=NULL;
}

static std::map <std::string, uint8_t> get_keyid={ {"KEYNAME", 0}, {"KEYID", 1}, {"KEYDESC",3} };
///The constructor to be used when keys are declared within data classes.
key::key(const std::string &line)
{
	std::cerr << line << '\n';
	std::vector <std::string> fields=split(line,' ');
	for (size_t x=0; x<fields.size(); ++x){
		std::vector <std::string> key_pair=split(fields[x],'=');
		uint8_t KEYID=get_keyid[key_pair[0]];
		std::cerr << "("<< fields[x] << ")=" << int(KEYID) << '\n';
		std::cerr << key_pair[0] << " and " << key_pair[1] << '\n';
/*		switch(KEYID){
			case 0:
				keyname_=new KEYNAME(key_pair[1].c_str() );
				break;
			case 1:
				keynum_=new KEYNUM(atoi(key_pair[1].c_str()));
				break;
			case 2:
				keydesc_=new KEYDESC(key_pair[1].c_str() );
				break;
			default:
				std::cerr << __FILE__ << ":" << __LINE__ << ": error: unrecognized KEYID.\n";
				exit(0);
		}*/
	}
	instance_=0;
	offset_=nokey;
//	size_of_f=&size_of_char7;
//	to_text_f=&to_text_char7;
//	from_tx_f=&from_tx_char7;
}
	
const char * key::get_name(void) const { return keyname_->value(); };
const uint8_t key::get_num(void) const { return *keynum_->value(); };
const std::string * key::get_desc(void) const { return keydesc_->value(); };
void * key::value(void) { return value_; };
void key::set_offset(const uint8_t &offset){
	offset_=offset;
};

key& key::operator=(key const& that){
	size_of_f=that.size_of_f;
	from_tx_f=that.from_tx_f;
	to_text_f=that.to_text_f;
	value_=that.value_;
	instance_=that.instance_;
	keyname_=that.keyname_;
	keynum_=that.keynum_;
	keydesc_=that.keydesc_;
	return *this;
}
