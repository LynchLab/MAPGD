#include "quartet.h"
/*
count_t major(const quartet_t);
count_t minor(const quartet_t);
count_t error1(const quartet_t);
count_t error2(const quartet_t);*/

count_t count(const quartet_t &c){
	return c.base[0]+c.base[1]+c.base[2]+c.base[3];
}
