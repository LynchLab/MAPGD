#include "base.h"

Base::Base (){
	base=5;
}

Base::Base (const char &c)
{
	base=Base::ctob(c);
}

Base::Base (const gt_t &g)
{
	base=g;
}

//These may be changed over to bit masks, in which case a lot of the 
//program will need to change.
gt_t Base::ctob(const char &c)
{
	switch(c){
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		case 'N':
			return 4;
		case 'a':
			return 0;
		case 'c':
			return 1;
		case 'g':
			return 2;
		case 't':
			return 3;
		case 'n':
			return 4;
		//R-V treat at of bits 5,6,7,8 
		case 'R':
			return 80;
		case 'Y':
			return 160;
		case 'K':
			return 192;
		case 'M':
			return 48;
		case 'W':
			return 144;
		case 'S':
			return 96;
		case 'B':
			return 224;
		case 'D':
			return 208;
		case 'H':
			return 176;
		case 'V':
			return 112;
		default :
			return 4;
			
	}
}

char Base::btoc(const gt_t &b) 
{
	switch(b){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
		case 4:
			return 'N';
		case 5:
			return 'a';
		case 6:
			return 'c';
		case 7:
			return 'g';
		case 8:
			return 't';
		case 9:
			return 'n';
		//R-V treat at of bits 5,6,7,8 
		case 80:
			return 'R';
		case 160:
			return 'Y';
		case 192:
			return 'K';
		case 48:
			return 'M';
		case 144:
			return 'W';
		case 96:
			return 'S';
		case 224:
			return 'B';
		case 208:
			return 'D';
		case 176:
			return 'H';
		case 112:
			return 'V';
		default :
			return 'N';
			
	}
	
}

std::istream& operator >> (std::istream& in, Base& x)
{
	char c; 
	in >> c;
	x.base=Base::ctob(c);
	return in;
}

std::ostream& operator<< (std::ostream& out, const Base& x)
{
	out << Base::btoc(x.base);
	return out;
}

/*
size_t Base::size(void) const 
{
	return sizeof(gt_t);
}*/
