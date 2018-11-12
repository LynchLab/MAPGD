/* synonym for population? */

#ifndef ALLELE_STAT_H
#define ALLELE_STAT_H

#include <iostream>
#include <iomanip>
#include <cfloat>
#include <sstream>
#include <vector>

#include "typedef.h"
#include "base.h"
#include "data.h"

///	Summary statistics from the allele command.
/**
 */
class Allele : public Indexed_data { 
private:
	void read(std::istream& str);
	///! the write function must be ? by the child class.
	void write(std::ostream& str) const;
	static const Registration registered;
	static Data * create(const std::vector <std::string> & Columns){
		return new Allele();
	};
public:
//	void write_binary (std::ostream& out) const;
//	void read_binary (std::istream& in);

	char delim;		//!< the delimiter used when reading/writing the class in text mode.	

	count_t excluded;	//!< A count of the number of samples that were excluded due to filtering criteria.
	bool pooled;		//!< Inferred from pooled or labeled sequencing?

	Allele();		//

	float_t freq;		//!< frequency of major allele.

	gt_t ref;		//!< identity of ref allele.
	gt_t minor;		//!< identity of minor allele.
	gt_t major;		//!< identity major allele.
	gt_t e1;		//!< identity of error1
	gt_t e2;		//!< identity of error2.

	float_t error;		//!< ml error rate.

	float_t null_error;	//!< error rate assuming Major allele monomorphism.
	float_t null_error2;	//!< error rate assuming Minor allele monomorphism.
	count_t coverage;	//!< population coverage.
	float_t efc; 		//!< number of 'effective' chromosomes in the sample.

	float_t ll, monoll, hwell; //!< log likelihoods.

	float_t MM; 		//!< frequency of genotype MM in the population.
	float_t Mm; 		//!< frequency of genotype Mm in the population.
	float_t mm; 		//!< frequency of genotype mm in the population.

	count_t N;		//!< number of individual at the site.
	float_t f;		//!< HW statistic.
	float_t h;		//!< heterozygosity.

	float_t bias;		//!< reference bias.
	float_t pbias;		//!< probability of reference bias.
	bool print_bias;	//!< probability of reference bias.

	float_t gof; 		//!< gof statistic.

	Allele & operator=(const Allele &);				//!< use the = operator to assign Allele.


	/// Allele needs to know whether to read the optional heterozygosity bias columns
	Allele(std::vector<std::string>);

	/*(0){}
	{
		std::cerr << "done!\n";
	};*/
	std::string header(void) const;	//!< returns first two lines of a file.
	static const std::string table_name;	//!< The destination table in the database.
	static const std::string file_name;	//!< Default file extension.
	static const bool binary;	        //!< Default file extension.
	size_t size(void) const;		//!< The size of the class in bytes.

	const std::string get_file_name(void) const;	//!< Returns the default file extension.
	const std::string get_table_name(void) const;	//!< Returns the destination table name.
	const bool get_binary(void) const;	        //!< Returns the destination table name.

	const std::string sql_header(void) const;	//!< Returns a string to create an SQL table.
	const std::string sql_column_names(void) const; //!< Returns the column names in the SQL table.
	const std::string sql_values(void) const;	//!< Returns a string to insert values in an SQL table.
};

#endif
