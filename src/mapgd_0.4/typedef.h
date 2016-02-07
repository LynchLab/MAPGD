#ifndef TYPEDEF_H_
#define TYPEDEF_H_

#include <cstddef>		//Included for ?
#include <stdint.h>		//Included for uint32_t, uint16_t
#include <cstdlib>		//Becuse why not?
#include <math.h>		//Included for float_t.

typedef uint16_t count_t;	//Should be used to specify deapth of coverage only.  
				//Once it is fixed it should be changed ot uint8_t or uint16_t.

typedef long double real_t;	//most ?

typedef uint16_t id0_t;		//specifies scaffold. 		limits to 65,536 scaffolds.

typedef uint32_t id1_t;		//specifies bp locations.	limits to 4,294,967,296 bp per scaffold.

typedef int64_t id1_off_t;	//specifies an offset		specifies a distance between bp, must be
				//				able to adress -id1_t to +id1_t.
	
typedef uint8_t gt_t;		//specifies a genotype.		limit 128.

#define CNT_MAX	65535
#define ID0_MAX	65535
#define ID1_MAX	4294967295
#define	GT_MAX	255
#define SQL_LINE_SIZE	1024
#endif
