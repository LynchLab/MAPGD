/*These are the commands that mapgd can run. 
 * Each one of them is a stand alone command line program the 
 * mapgd calls in order to something or other... Perserve the 
 * name space. Save the children. I don't know. 
*/



/** \defgroup Command Mapgd commands
 *  These are the basic commands that can be called by the main mapgd 
 *  program.
 */
 
#ifndef _COMMANDS_H_
#define _COMMANDS_H_

#include "allele_cmd.h"
#include "estimate_pooled.h"
#include "allele.h"
#include "estimate_fst.h"
#include "relatedness.h"
#include "proview.h"
#include "filter.h"
#include "fastview.h"
//#include "meta_filter.h"
#include "filter_pool.h"
#include "filter_genotype.h"
//#include "filter_pro.h"
#include "map2genotype.h"
#include "linkage_disequilibrium.h"
#include "sam2idx.h"
#include "make_vcf.h"
#include "make_vcf2.h"
#include "read_vcf.h"
#include "mapgd_help.h"
#include "relatedness_test.h"
#include "test_keys.h"
#include "read_pheno.h"
#include "read_bed.h"
#include "surface.h"

//#include "simulate.h"

#ifndef NOSQL
#include "writesql.h"
#include "readsql.h"
#endif 

#endif
