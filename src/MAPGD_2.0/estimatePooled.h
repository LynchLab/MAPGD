#ifndef _ESTIMATOR_HPP_
#define _ESTIMATOR_HPP_	

#include	<stdio.h>
#include 	<math.h>
#include 	<sys/types.h>
#include	<stdlib.h>
#include 	<iostream>
#include 	"interface.h" 
#include 	"proFile.h" 

int estimatePooled(int, char **);

typedef struct allele_stat { 
	float_t freqmajor;
	float_t freqminor;
	float_t errorrate;
	float_t coverage;
	float_t ll; 

} allele_stat_t;

allele_stat_t estimate(quartet_t);
#endif 
