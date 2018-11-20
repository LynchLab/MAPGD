/* A small set of functions so I can pretend I'm writing in python. */

#include "stream_tools.h"

#ifdef POSIX

bool 
check_stream(std::istream *s)
{
	return false;
	return true;
}
#endif


//This function is based off of code written by user Drise on https://stackoverflow.com/
/*
istream &
parse_on_fail (istream &in, double &x, bool neg) {
    const char *exp[] = { "", "inf", "NaN" };
    const char *e = exp[0];
    int l = 0;
    char inf[4];
    char *c = inf;
    if (neg) *c++ = '-';
    in.clear();
    if (!(in >> *c).good()) return *this;
    switch (*c) {
    case 'i': e = exp[l=1]; break;
    case 'N': e = exp[l=2]; break;
    }
    while (*c == *e) {
        if ((e-exp[l]) == 2) break;
        ++e; if (!(in >> *++c).good()) break;
    }
    if (in.good() && *c == *e) {
        switch (l) {
        case 1: x = std::numeric_limits<double>::infinity(); break;
        case 2: x = std::numeric_limits<double>::quiet_NaN(); break;
        }
        if (neg) x = -x;
        return *this;
    } else if (!in.good()) {
        if (!in.fail()) return *this;
        in.clear(); --c;
    }
    do { in.putback(*c); } while (c-- != inf);
    in.setstate(std::ios_base::failbit);
    return *this;
}*/
