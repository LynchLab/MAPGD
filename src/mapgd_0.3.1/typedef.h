#ifndef TYPEDEF_H_
#define TYPEDEF_H_

#define VERSION	"0.4"		//Since binary compatibility will likely depend on having the sizes defined in this file stay
				//constant, the version number needs to be incremented whenever anything in this file changes.
				//The program should always be able to open files after VERSION 1.0

#include <cstddef>		//Included for ?
#include <stdint.h>		//Included for uint32_t, uint16_t
#include <cstdlib>		//Becuse why not?
#include <math.h>		//Included for float_t.

typedef uint8_t count_t;	//Should be used to specify deapth of coverage only.  
				//Once it is fixed it should be changed ot uint8_t or uint16_t.

typedef long double real_t;	//most ?

typedef uint16_t id0_t;		//specifies scaffold. 		limits to 65,536 scaffolds.

typedef uint32_t id1_t;		//specifies bp locations.	limits to 4,294,967,296 bp per scaffold.

typedef int64_t id1_off_t;	//specifies an offset		specifies a distance between bp, must be
				//				able to adress -id1_t to +id1_t.
	
typedef uint8_t gt_t;		//specifies a genotype.		limit 128.

#endif
