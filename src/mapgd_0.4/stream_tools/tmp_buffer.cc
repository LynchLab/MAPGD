#include "tmp_buffer.h"

Tmp_streambuf::Tmp_streambuf()
{
	_opened=false;
	buffered=false;
}

Tmp_streambuf *
Tmp_streambuf::open(std::istream **dest, std::istream *src)
{
//	std::cerr << "Opening!\n";
	if ( (*dest)->good() && src->good() )
	{
//		std::cerr << "Good!\n";
		_mode = READ;
		_opened = true;
		_src=src;
		_dest=dest;
		buffered=false;
		reread=false;
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ << "Could not open i/o streams." << std::endl;
	}
	return this;
}

int 
Tmp_streambuf::underflow() 
{
//	std::cerr << "Underflow!\n";
	if (_opened) {
		if(!reread){
			c=_src->rdbuf()->sbumpc();
			if (buffered) _buffer << c;
//			std::cerr << buffered\n";
		} else {
			_buffer.get(c);
			if (!_buffer.good() ){
				c=_src->rdbuf()->sbumpc();
				*_dest=_src;
			}
		}
		//std::cerr << c << std::endl;
		setg(&c, &c, &c+1);
		return c;
	}
	else return EOF;
}

void 
Tmp_streambuf::clear_read(void)
{
	//std::cerr << "Clear!\n";
	_buffer.str(std::string());
	//std::cerr << "Good? " << _buffer.good() << std::endl;
	//_
}

void Tmp_buffer::clear_read(void)
{
	_spool.clear_read();
}

bool 
Tmp_streambuf::is_open(void)
{
	return _opened;
}

Tmp_streambuf *
Tmp_streambuf::close(void)
{
//	std::cerr << "Close!\n";
//	_opened=false;
	return (Tmp_streambuf*)0;
}

Tmp_streambuf::~Tmp_streambuf()
{
//	std::cerr << "Something is happening!\n";
	close(); 
}

Tmp_buffer::Tmp_buffer (std::istream **dest, std::istream *src) : std::istream(NULL) 
{
	*dest=this;
	init(&_spool);
	_spool.open(dest, src);
}

void
Tmp_buffer::open (std::istream **dest, std::istream *src) 
{

	*dest=this;
	init(&_spool);
	_spool.open(dest, src);
}
