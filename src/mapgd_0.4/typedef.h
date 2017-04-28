#ifndef _TYPEDEF_H_
#define _TYPEDEF_H_

#include <cstddef>		//Included for ?
#include <stdint.h>		//Included for uint32_t, uint16_t
#include <cstdlib>		//Becuse why not?
#include <math.h>		//Included for float_t.
#include <float.h>		//For the limits of floating types.

/* Some people might consider these types of typedefs "cutesy", in that they 
 * don't do much, but we really want to limit the size of some of these objects
 * in memory, and so we want to find the smallest data ...
 */
typedef uint16_t count_t;	//Should be used to specify depth of coverage only.  
				//Once it is fixed it should be changed to uint8_t or uint16_t.

typedef double real_t;		//most ?
typedef uint16_t id0_t;		//specifies scaffold. 		limits to 65,536 scaffolds.
typedef uint32_t id1_t;		//specifies bp locations.	limits to 4,294,967,296 bp per scaffold.
typedef int64_t id1_off_t;	//specifies an offset		specifies a distance between bp, must be
typedef uint8_t gt_t;		//specifies a genotype.		limit 128.

#define CNT_MAX	65535
#define ID0_MAX	65535
#define ID1_MAX	4294967295
#define REAL_MAX LDBL_MAX
#define	GT_MAX	255
#define SQL_LINE_SIZE	1024


#endif
