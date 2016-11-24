#ifndef _READ_VCF_H_
#define _READ_VCF_H_	

#include "allele.h"
#include "locus.h"
#include "map-file.h"
#include "interface.h"
#include "vcf-file.h"

#ifdef MPI
#include <ciso646>
#include <mpi.h>
#endif

int read_vcf(int, char **);

#endif 
