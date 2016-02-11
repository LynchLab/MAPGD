#ifndef _LINKAGE_H_
#define _LINKAGE_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#include "../typedef.h"
#include "data.h"
#include "../stream-tools.h"

///	A class to store linkage information. 
/** Many of these data classes need to have things moved over to static 
  * functions. The problem popped up because I just copied over all the output 
  * people wanted into values in the class, but really a lot of it can be 
  * determined from ... bah. I'm getting bored of talking to myself. Someday
  * if I ever get any help I may actually make these comments useful. Until 
  * please enjoy this link to a 
<a href="https://www.youtube.com/watch?v=VXa9tXcMhXQ">Kraftwork video</a> to see if Doxygen automatically 
  * incorporates HTML tags right. Listening to Kraftwerk improves your coding 
  * skillz. It's totally a fact.
  * I'll start fixing allele stat now.  
 */
class Linkage : public virtual Data { 
private:
	void read(std::istream& str);
	///! the write function must be ? by the child class.
	void write(std::ostream& str) const;
	static const Registration registered;
	static Data * create(const std::vector <std::string> &columns){
		return new Linkage(columns);
	};
	//We are comparing site X to site Y, so we need a total of four ids.
	id0_t id0_X_;	//!< The scaffold of site X
	id1_t id1_X_;	//!< The pos of site X

	id0_t id0_Y_;	//!< the scaffold of site Y 
	id1_t id1_Y_;	//!< the pos of site Y

	float_t p_;		//!< freq_major site numero uno
	float_t q_;		//!< freq_major site numero dos

	float_t D_;		//!< the magnitude of the ld, j0 
	float_t Ni_;		//!< The number of individuals used in the calculation
	float_t fit_;		//!< fit statistic (thar be ld here!)
	float_t null_;		//!< null statistic (i.e. no ld)
public:
	char delim;		//!< the delimiter used when reading/writing the class in text mode.	

	Linkage();	
	Linkage(const std::vector <std::string> &) : Linkage(){};	//!< delegating the constructor ftw.	

	static const std::string file_name;				//!< The default extension for files.
	static const std::string table_name;				//!< Destination table in Db.

	const std::string header(void) const;
	size_t size(void) const;

	const std::string get_file_name(void) const;				//!< The default extension for files.
	const std::string get_table_name(void) const;				//!< The default extension for files.
	const std::string sql_header(void) const;				//!< Destination table in Db.
	const std::string sql_column_names(void) const;				//!< Destination table in Db.
	const std::string sql_values(void) const;				//!< Destination table in Db.

	void set_D(const float_t &);
	void set_fit(const float_t &);
	void set_null(const float_t &);
	void set_Ni(const float_t &);
	void set_p(const float_t &);
	void set_q(const float_t &);

	float_t ll (void) const;
	float_t dist (void) const;
	float_t D (void) const;
	float_t Dmax (void) const;
	float_t Dmin (void) const;
	float_t Dsq (void) const;
	float_t Dprime (void) const;
	float_t r (void) const;
	float_t rsq (void) const;
	float_t adj_D (void) const;
	float_t adj_Dsq (void) const;
	float_t adj_Dprime (void) const;
	float_t adj_r (void) const;
	float_t adj_rsq (void) const;
	float_t adj_Dmax (void) const;
	float_t adj_Dmin (void) const;
};

#endif
