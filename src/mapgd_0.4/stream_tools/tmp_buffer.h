#include <iostream>
#include <sstream>
#include <string>
#include <cstring>

#ifndef TMP_BUFFER_H
#define TMP_BUFFER_H

#ifndef EOF
#define EOF  -1
#endif

#ifndef READ
#define READ std::ios::in
#endif

class Tmp_streambuf : public std::streambuf {
	private:
	bool _opened;
	int _mode;
	char c;
	std::istream *_src;
	std::istream **_dest;
	std::stringstream _buffer;
	public:
	bool buffered, reread;
	Tmp_streambuf();
	Tmp_streambuf * open(std::istream **dest, std::istream *src);
	virtual int underflow();
	bool is_open(void);
	Tmp_streambuf * close(void);
	~Tmp_streambuf();
	void clear_read();
};

/***/
class Tmp_buffer : public std::istream {

	private:

	Tmp_streambuf _spool;

	protected:
	public:
	
	Tmp_buffer (void) : std::istream(NULL){};
	~Tmp_buffer (void){};

	Tmp_buffer (std::istream **dest, std::istream *src);
	void open(std::istream **dest, std::istream *src);
	void buffer_on (void) {_spool.buffered=true;};
	void buffer_off (void) {_spool.buffered=false;};
	void reread_on (void) {_spool.reread=true;};
	void reread_off (void) {_spool.reread=false;};
	void clear_read();
};
#endif
