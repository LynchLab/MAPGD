#ifndef _READ_BED_H_
#define _READ_BED_H_	

#include "allele.h"
#include "locus.h"
#include "map_file.h"
#include "interface.h"
#include "file_index.h"
#include "bed_file.h"
#include "system.h"

#ifndef NOLZ4
#include "state.h"
#endif 

int read_bed(int, char **);

#endif 
