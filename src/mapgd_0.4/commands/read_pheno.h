#ifndef _READ_PHENO_H_
#define _READ_PHENO_H_	

#include "phenotype.h"
#include "map_file.h"
#include "interface.h"
#include "plink_pheno.h"

#ifdef MPI
#include <ciso646>
#include <mpi.h>
#endif

int read_pheno(int, char **);

#endif 
