#ifndef _MAKE_VCF_H_
#define _MAKE_VCF_H_	

#include "allele.h"
#include "locus.h"
#include "map_file.h"
#include "interface.h"
#include "vcf_file.h"

#ifdef MPI
#include <ciso646>
#include <mpi.h>
#endif

int make_vcf(int, char **);

#endif 
