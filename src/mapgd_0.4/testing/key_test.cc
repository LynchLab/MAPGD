#include "../datatypes/row.h"
#include "../typedef.h"
#include "../datatypes/base_keys.h"

#include <math.h>
#include <cstdio>

int main (int argc, char *argv[]){
	std::cerr << "row this_row\n";
	row this_row(std::cin);
	std::cerr << "get_row\n";
	int x=10;
	void * ptr=&x;
	int& y=*(int *)(ptr);
	if ( not( get_row(std::cin, this_row) ) ){
		std::cerr << "i stream not open\n";
		exit(0);
	}
	std::cerr << "real_key\n";
	real_key freq_key;
	std::cerr << "freq_key=\n";
	freq_key=this_row.get_key(keyid::freq);
	std::cerr << "freq_key.value\n";
	real_t *X=freq_key.value();
	real_t lnsum=0;
	while (get_row(std::cin, this_row) ){
		lnsum+=log(*X);
	}
	return 0;
}
