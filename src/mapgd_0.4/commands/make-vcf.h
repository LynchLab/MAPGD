#ifndef _VCF_H_
#define _VCF_H_	

#include "allele.h"
#include "locus.h"
#include "map-file.h"
#include "interface.h"

#ifdef MPI
#include <ciso646>
#include <mpi.h>
#endif

int vcf(int, char **);

#endif 
