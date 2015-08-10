#ifndef _ALLELE_H_
#define _ALLELE_H_

#include "key.h"
#include "../typedef.h"

///	A class to store population specific information. May be moved over to population.
typedef struct allele {
	real_t frequency[3], error;
	gt_t base[4];
	allele & operator=(const allele &);	//!< use the = operator to assign allele.
	
} __attribute__((packed)) allele_t ;

inline real_t get_MM(const allele_t &a){ return a.frequency[0];} //!< frequency of genotype MM in the popualtion.
inline real_t get_Mm(const allele_t &a){ return a.frequency[1];} //!< frequency of genotype Mm in the popualtion.
inline real_t get_mm(const allele_t &a){ return a.frequency[2];} //!< frequency of genotype mm in the popualtion.

inline gt_t get_major(const allele_t &a) { return a.base[0]; }	//!< idenity of major allele.
inline gt_t get_minor(const allele_t &a) { return a.base[1]; }	//!< idenity of minor allele.
inline gt_t get_error1(const allele_t &a) { return a.base[2]; }	//!< idenity of an error allele.
inline gt_t get_error2(const allele_t &a) { return a.base[3]; }	//!< idenity of an error allele.

inline real_t get_freq(const allele_t &a) {return a.frequency[0]+a.frequency[1]/2.;}; //!< returns the frequency of the major allele
inline real_t get_f(const allele_t &a) {
	real_t p=get_freq(a);
	return 1-a.frequency[1]/(2.*p*(1.-p) );
} //!< returns a measure of departure from Hardy-Weinberg equilibrium (the f statistic)

#endif
/*
	real_t &MM=frequency[0];
	real_t &Mm=frequency[1];
	real_t &mm=frequency[2];

	gt_t &major=base[0];
	gt_t &minor=base[1];
	gt_t &error1=base[2];
	gt_t &error2=base[3];*/
