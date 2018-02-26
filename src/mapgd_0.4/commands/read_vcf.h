#ifndef _READ_VCF_H_
#define _READ_VCF_H_	

#include "allele.h"
#include "locus.h"
#include "map_file.h"
#include "interface.h"
#include "file_index.h"
#include "vcf_file.h"

#ifndef NOLZ4
#include "state.h"
#endif 

#ifdef MPI
#include <ciso646>
#include <mpi.h>
#endif

int read_vcf(int, char **);

#endif 
