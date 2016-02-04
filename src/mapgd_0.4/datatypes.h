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
#include "allele_stat.h"	//FORMALLY allele_stat.
#include "pooled_data.h"	//FORMALLY allele_stat.
#include "data_types/sample_name.h"	//FORMALLY allele_stat.
#include "genotype.h"

#endif 
