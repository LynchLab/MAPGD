/* An bcf  */

#ifndef _BCF_FILE_H_
#define _BCF_FILE_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#include "typedef.h"
#include "datatypes.h"
#include "raw.h"
#include "stream_tools.h"

#include "map_file.h"
#include "locus.h"
#include "bcf2pro.h"

/// Because of the god awful mess that are vcf header lines.
/** This is likely to become some form of container to handle moving data into and out of rows of map file.
 */

class Bcf2pro_file : public Indexed_file <Locus> {
	using Indexed_file<Locus>::file_index_;

	void read_text(Locus &);	//!< Read file in text mode.			DONE
	void read_text_profile(Locus &locus){Indexed_file<Locus>::read_text(locus);};	//!< Read file in text mode.			DONE

	bool profile_;

	int columns_, offset_;

	using Indexed_file<Locus>::write_text;	//!< Write in text mode.	DONE
	using Indexed_file<Locus>::out_;	//(const std::ios_base::openmode &);
	using Indexed_file<Locus>::in_;	//(const std::ios_base::openmode &);

	using Base_file::open_;
	using Base_file::table_open_;
public:
	Bcf2pro_file(){profile_=false;};
	Bcf2pro_file(bool call){profile_=call;};
	void set_mpileup(const int &, const int &);
/*	Mpileup_file(){
		open_=false;
		table_open_=false;
		read_=false;
		write_=false;
		delim_column_='\t';
		binary_=false;
		filename_="";
		concatenated_=false;
		indexed_=false;
	};			//!< default constructor*/
	using Indexed_file<Locus>::open;
	using Indexed_file<Locus>::set_index;	//!< Sets the File_index.		
	using Indexed_file<Locus>::get_index;	//!< Returns the File_index.

	/*! \brief Returns the position in the file.
	 *
	 * Positions are guarnteed to be unique and increasing for each row. 
	 */

	using Indexed_file<Locus>::get_pos;
	using Indexed_file<Locus>::write_header;	//!< Writes a file header.
	
	Locus read_header_profile(void){return Indexed_file<Locus>::read_header();}	//!< Writes a file header.
	Locus read_header(void);			//!< Reads a file header.
};

#endif
