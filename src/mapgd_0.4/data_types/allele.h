/* synonym for population? */

#ifndef ALLELE_STAT_H
#define ALLELE_STAT_H

#include <iostream>
#include <iomanip>
#include <cfloat>
#include <sstream>
#include <vector>

#include "../typedef.h"
#include "../base.h"
#include "data.h"

///	A class to store population specific information. May be moved over to population.
/** This is likely to become some form of container to handel moving data into and out of rows of map file.
 */
class Allele : public Data { 
private:
	void read(std::istream& str);
	///! the write function must be ? by the child class.
	void write(std::ostream& str) const;
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Allele();
	};
public:
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	

	id0_t id0;		//!< the scaffold identifer of the allele.
	id1_t id1;		//!< the bp location of the allele.

	count_t excluded;	//!< A count of the number of samples that were excluded due to filtering criteria.
	bool pooled;		//!< Infered from pooled or labeled sequencing?

	Allele();		//

	float_t freq;		//!< frequency of major allele.

	gt_t ref;		//!< idenity of ref allele.
	gt_t minor;		//!< idenity of minor allele.
	gt_t major;		//!< idenity major allele.
	gt_t e1;		//!< idenity of error1
	gt_t e2;		//!< idenity of error2.

	float_t error;		//!< ml error rate.

	float_t null_error;	//!< error rate assuming monomorphism.
	count_t coverage;	//!< population coverage.
	float_t efc; 		//!< number of 'effective' chromosomes in the sample.

	float_t ll, monoll, hwell; //!< loglikelihoods.

	float_t MM; 		//!< frequency of genotype MM in the population.
	float_t Mm; 		//!< frequency of genotype Mm in the population.
	float_t mm; 		//!< frequency of genotype mm in the popualtion.

	count_t N;		//!< number of individual at the site.
	float_t f;		//!< HW statistic.
	float_t h;		//!< heterozygosity.

	float_t gof; 		//!< gof statistic.

	Allele & operator=(const Allele &);				//!< use the = operator to assign Allele.

	/**
	 * \defgroup basic_data
	 * @{
	 */

	///! Required constructor 
	Allele(std::vector<std::string>) : Allele(){};
	const std::string header(void) const;		//!< names of columns.
	static const std::string table_name;	//!< The destination table in the Db.
	static const std::string file_name;	//!< Default file extention.
	size_t size(void) const;		//!< The size of the ? in bytes.

	const std::string get_file_name(void) const;				//!< The dafualt extention for files.
	const std::string get_table_name(void) const;				//!< The dafualt extention for files.
	const std::string sql_header() const;
	const std::string sql_column_names(void) const;
	const std::string sql_values(void) const;
	 // @}
};

#endif
