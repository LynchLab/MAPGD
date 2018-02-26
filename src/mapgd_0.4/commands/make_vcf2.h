#ifndef _MAKE_VCF2_H_
#define _MAKE_VCF2_H_	

#include "allele.h"
#include "locus.h"
#include "map_file.h"
#include "interface.h"
#include "vcf_file.h"

#ifndef NOLZ4
#include "state.h"
#endif

#ifdef MPI
#include <ciso646>
#include <mpi.h>
#endif

int make_vcf2(int, char **);

#endif 
