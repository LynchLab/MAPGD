/*These are the commands that mapgd can run. 
 * Each one of them is a stand alone command line program the 
 * mapgd calls in order to something or other... Perserve the 
 * name space. Save the children. I don't know. 
*/



/** @defgroup Command Mapgd commands
 *  These are the basic commands that can be called by the main mapgd 
 *  program.
 */
 
#ifndef _COMMANDS_H_
#define _COMMANDS_H_

//#include "compare-pooled.h"
#include "estimate-pooled.h"
#include "estimate-individual.h"
#include "estimate-fst.h"
#include "relatedness.h"
#include "proview.h"
//#include "convert.h"
#include "filter.h"
#include "map2genotype.h"
//#include "PopLD.h"
#include "sam2idx.h"
//#include "sql/readsql.h"
#include "sql/writesql.h"

#endif
