/* synonym for population? */

#ifndef ALLELE_STAT_H
#define ALLELE_STAT_H

#include <iostream>
#include "typedef.h"
#include <cfloat>
#include <iomanip>

///	A class to store population specific information. May be moved over to population.
/** This is likely to become some form of container to handel moving data into and out of rows of map file.
 */
class allele_stat { 
private:
public:

	char delim;		//!< the delimiter used when reading/writing the class in text mode.	

	std::string id0;	//!< the scaffold identifer of the allele.
	id1_t id1;		//!< the bp location of the allele.

	size_t excluded;	//!< A count of the number of samples that were excluded due to filtering criteria.
	bool pooled;		//!< Infered from pooled or labeled sequencing?

	allele_stat();		//

	map_double_t freq;		//!< frequency of major allele.

	gt_t minor;		//!< idenity of minor allele.
	gt_t major;		//!< idenity major allele.
	gt_t e1;		//!< idenity of minor allele.
	gt_t e2;		//!< idenity major allele.

	map_double_t error;		//!< ml error rate.

	map_double_t null_error;	//!< error rate assuming monomorphism.
	map_double_t coverage;	//!< population coverage.
	map_double_t efc; 		//!< number of 'effective' chromosomes in the sample.

	map_double_t ll, monoll, hwell; //!< loglikelihoods.

	//FOR LABELED SEQUENCING ONLY!
	map_double_t MM; 		//!< frequency of genotype MM in the population.
	map_double_t Mm; 		//!< frequency of genotype Mm in the population.
	map_double_t mm; 		//!< frequency of genotype mm in the popualtion.

	map_double_t N;		//!< number of individual at the site.
	map_double_t f;		//!< HW statistic.
	map_double_t h;		//!< heterozygosity.

	map_double_t gof; 		//!< gof statistic.

	allele_stat & operator=(const allele_stat &);				//!< use the = operator to assign allele_stat.
	friend std::ostream& operator<< (std::ostream&, const allele_stat&);	//!< use the << operator to write allele_stat.
};

#endif
