/* An vcf  */

#ifndef _VCF_FILE_H_
#define _VCF_FILE_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#ifndef NOHTS
#include <htslib/hts.h>
#include <htslib/vcf.h>

#ifndef NOLZ4
#include "state.h"
#endif

#include "external-file.h"
#include "external-data.h"

#include "typedef.h"
#include "datatypes.h"
#include "raw.h"
#include "stream_tools.h"

#ifndef bcf_get_format_int32
#define bcf_get_format_int32(hdr,line,tag,dst,ndst)	bcf_get_format_values(hdr,line,tag,(void **)(dst),ndst,BCF_HT_INT)
#endif

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

//	bcf_init
//	bcf_destroy
//	bcf_read(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)
//	bcf_write(htsFile *fp, const bcf_hdr_t *h, bcf1_t *v)

	//! number of individuals sapmled

public:
	Vcf_data ();

	//These should all be private, make vcf_file freind?
	bcf1_t *record_;
	bcf_hdr_t *header_;
	size_t sample_size_;

	void set_header(const File_index &, const std::vector <std::string> &);

	void put (const Data *, ...);
	void get (Data *, ...) const;

	void put (const File_index &, const Allele &, const Population &);
	id1_t get (const File_index &, Population &) const;

#ifndef NOLZ4
	id1_t get (State &) const;
#endif 
	std::vector<std::string> get_sample_names (void) const;
	File_index get_index (void) const;

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
	htsFile *file_;
public:
	void open(const std::ios_base::openmode &);
	void open(const char *, const std::ios_base::openmode &);
	void close(void);
	void read (Vcf_data &);
	void write (const Vcf_data &);
	void write_header (const Vcf_data &);
	Vcf_data read_header();
};

/*
typedef struct {
    int32_t n[3];           // n:the size of the dictionary block in use, (allocated size, m, is below to preserve ABI) 
    bcf_idpair_t *id[3];
    void *dict[3];          // ID dictionary, contig dict and sample dict
    char **samples;	    // Presumably sample names.
    bcf_hrec_t **hrec;      // Jesus fucking crist.
    int nhrec, dirty;       // 
    int ntransl, *transl[2];    // for bcf_translate()
    int nsamples_ori;           // for bcf_hdr_set_samples()
    uint8_t *keep_samples;
    kstring_t mem;
    int32_t m[3];          // m: allocated size of the dictionary block in use (see n above)
} bcf_hdr_t;
*/

#endif
#endif
