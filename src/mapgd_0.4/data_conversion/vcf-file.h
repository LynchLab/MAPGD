/* An vcf  */

#ifndef _VCF_FILE_H_
#define _VCF_FILE_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#include "htslib/vcf.h"

#include "external-file.h"
#include "external-data.h"

#include "typedef.h"
#include "datatypes.h"
#include "raw.h"
#include "stream-tools.h"

/// Because of the god awful mess that are vcf header lines.
/** This is likely to become some form of container to handle moving data into and out of rows of map file.
 */

class Info {
public:
	enum Key {
		AA,
		BQ,
		CIGAR,
		DB,
		DP,
		END,
		H2,
		H3,
		MQ,
		MQ0,
		NS,
		SB,
		SOMATIC,
		VALIDATED,
		G
	};

};

class Vcf_data : public External_data {
private:
	bcf1_t *record;
	bcf_hdr_t *header;
//	bcf_init
//	bcf_destroy
//	bcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
//	bcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)

public:
	void read (Allele& );
	void write (Allele& ) const;
	Vcf_data ();
	//The mandatory fields.
	std::string id;
	Base ref;
	std::vector <std::string> alt;
	float_t qual;
	bool filter;
	std::vector <std::string> info;

	//! Ancestral allele
        Base AA;

	//! Allele count in genotypes, for each ALT allele, in the same order as listed
	std::vector <count_t> AC;

	//! Allele frequency for each ALT allele in the same order as listed: 
	//  use this when estimated from primary data, not called genotypes
	std::vector <float_t> AF;

 	//! Total number of alleles in called genotypes
        count_t AN;

	//! RMS base quality at this position
        float_t BQ;

	//! Cigar string describing how to align an alternate allele to the reference allele
	//cigar CIGAR;

 	//! dbSNP membership
	bool DB;

 	//! Combined depth across samples, e.g. DP=154
        count_t DP;

 	//! End position of the variant described in this record (esp. for CNVs)
	id1_t END;

 	//! Membership in hapmap2
	bool H2;

 	//! RMS mapping quality, e.g. MQ=52
	float_t MQ;

 	//! Number of MAPQ == 0 reads covering this record
	count_t MQ0;

        //! Number of samples with data
	count_t NS;

	//! Strand bias at this position
	float_t SB;

	//! Indicates that the record is a somatic mutation, for cancer genomics
        bool SOMATIC;

	//! Validated by follow-up experiment
        bool VALIDATED; 
};


class Vcf_file : public External_file <Vcf_data> {
private:
	using Base_file::in_;
	using Base_file::out_;
	using Base_file::open_;
	using Base_file::open_no_extention;

public:
	void open(const std::ios_base::openmode &);
	void open(const char *, const std::ios_base::openmode &);
	void read (Allele &allele);
	void write (const Allele &out) const;
};

#endif
