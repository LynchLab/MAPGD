/* This file should include all data types that are supported for 
 * reading and writing to files. Indexed files need to have the id0 
 * and id1 members (which specify their positions in the reference)
 * and all data types need a table_name, file_name, header(), >>, 
 * <<, size(), and a constructor from column names. See map-file.cc
 */

#ifndef DATA_TYPES_H_
#define DATA_TYPES_H_

#include "quartet.h"
#include "locus.h"
#include "data_types/allele.h"	//FORMALLY Allele.
#include "pooled_data.h"	//FORMALLY Allele.
#include "data_types/sample_name.h"	//FORMALLY Allele.
#include "genotype.h"
#include "data_types/linkage_data.h"

#endif 
