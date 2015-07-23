/* synonym for population? */

#ifndef ALLELE_STAT_H
#define ALLELE_STAT_H

#include <iostream>
#include "typedef.h"
#include <cfloat>
#include <iomanip>

class allele_stat { 
private:
public:

	char delim;		
//
	std::string id0;	//
	count_t id1;		//
//
	count_t excluded;	//
	bool pooled;		//Infered from pooled or labeled sequencing?

	allele_stat();		//

	float_t freq;		//frequency of major allele.

	count_t minor;		//idenity of minor allele.
	count_t major;		//idenity major allele.
	count_t e1;		//idenity of minor allele.
	count_t e2;		//idenity major allele.

	float_t error;		//ml error rate.

	float_t null_error;	//ml error rate.
	float_t coverage;	//population coverage.
	float_t efc; 		//number of 'effective' chromosomes in the sample.

	float_t ll, monoll, hwell; //loglikelihoods.

	//FOR LABELED SEQUENCING ONLY!
	float_t MM; 		//frequency of genotype MM in the population.
	float_t Mm; 		//frequency of genotype Mm in the population.
	float_t mm; 		//frequency of genotype mm in the popualtion.

	float_t N;		//number of individual at the site.
	float_t f;		//HW statistic.
	float_t h;		//heterozygosity.

	float_t gof; 		//gof statistic.

	allele_stat & operator=(const allele_stat &);
	friend std::ostream& operator<< (std::ostream&, const allele_stat&);
};

#endif
