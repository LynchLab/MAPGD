/* synonym for population? */

#ifndef RELATEDNESS_DATA_H_
#define RELATEDNESS_DATA_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>
#include <sstream>

#include "typedef.h"

///	A class to store population specific information. May be moved over to population.
/** This is likely to become some form of container to handel moving data into and out of rows of map file.
 */
class Relatedness { 
private:
public:
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	

	std::string X_;
	std::string Y_;

	id1_t sites;

	float_t e_X_[8], e_X_ll;
	float_t e_Y_[8], e_Y_ll;
	float_t f_X_, f_X_ll;
	float_t f_Y_, f_Y_ll;
	float_t theta_XY_, theta_XY_ll;
	float_t gamma_XY_, gamma_XY_ll;
	float_t gamma_YX_, gamma_YX_ll;
	float_t delta_XY_, delta_XY_ll;
	float_t Delta_XY_, Delta_XY_ll;
	float_t null_ll_, max_ll_;

	Relatedness();	
	Relatedness(const std::vector <std::string> &) : Relatedness(){};	
	Relatedness(const std::string &, const std::string &);	

	std::string header(void) const;
	size_t size(void) const;

	void set_X_name(const std::string &);
	void set_Y_name(const std::string &);

	friend std::ostream& operator<< (std::ostream&, const Relatedness&);	//!< use the << operator to write allele_stat.
	friend std::istream& operator>> (std::istream&, Relatedness&);		//!< use the >> operator to read allele_stat.
	static const std::string file_name;					//!< The dafualt extention for files.
	static const std::string table_name;					//!< Destination table in Db.
};

#endif
