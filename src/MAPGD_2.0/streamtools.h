#ifndef STREAMTOOLS_H_
#define STREAMTOOLS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

/*another tutorail code snipit*/
class opipe
{
public:
    opipe(std::ostream & dst, std::ostream & src)
        : src(src), sbuf(src.rdbuf(dst.rdbuf())) {}
//    ~opipe2() { src.rdbuf(sbuf); }
private:
    std::ostream & src;
    std::streambuf * const sbuf;
};

class ipipe
{
public:
    ipipe(std::istream & dst, std::istream & src)
        : src(src), sbuf(src.rdbuf(dst.rdbuf())) {}
//    ~ipipe2() { src.rdbuf(sbuf); }
private:
    std::istream & src;
    std::streambuf * const sbuf;
};

/*Thanks Viet for this code!*/
struct opiped {
    opiped(std::streambuf * buf, std::ostream & os)
    :os(os), old_buf(os.rdbuf(buf)) { }
    ~opiped() { os.rdbuf(old_buf); }

    std::ostream& os;
    std::streambuf * old_buf;
};

struct ipiped {
    ipiped(std::streambuf * buf, std::istream & is)
    :is(is), old_buf(is.rdbuf(buf)) { }
    ~ipiped() { is.rdbuf(old_buf); }

    std::istream& is;
    std::streambuf * old_buf;
};

std::vector<std::string> split(std::string s, char delim);
#endif
