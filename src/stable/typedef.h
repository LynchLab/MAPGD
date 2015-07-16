#ifndef TYPEDEF_H_
#define TYPEDEF_H_
#include <stdint.h>

#define VERSION	"1.0"

typedef uint32_t count_t;

#if __SIZEOF_FLOAT__ == 4
	typedef long double float_t;
#else
	typedef float float_t;
#endif 

#endif
