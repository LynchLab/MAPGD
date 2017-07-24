/* An bcf  */

#ifndef _BCF2PRO_H_
#define _BCF2PRO_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#include "typedef.h"
#include "datatypes.h"
#include "raw.h"
#include "stream_tools.h"

#include "map_file.h"

/// Because of the god awful mess that are vcf header lines.
/** This is likely to become some form of container to handle moving data into and out of rows of map file.
 */

class Bcf2pro : public Locus  {
private:

//	bcf_init
//	bcf_destroy
//	bcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
//	bcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)

public:
	Bcf2pro (){
//		Locus();
	};
	Bcf2pro (const std::vector <std::string> &str){
//..		Locus(str);
	};
	Bcf2pro (const count_t &){
//..		Locus(str);
	};


/*
	void set_header(const File_index &, const std::vector <std::string> &);

	void put (const Data *, ...);
	void get (Data *, ...) const;

	void put (const File_index &, const Locus &);
	void get (Locus &) const;

	//The mandatory fields.
	std::string id;
	Base ref;
	std::vector <std::string> alt;
	float_t qual;
	bool filter;
	std::vector <std::string> info;
*/
};

#endif
