/* synonym for population? */

#ifndef _VCF_H_
#define _VCF_H_

#include <iostream>
#include <cfloat>
#include <iomanip>
#include <vector>

#include "typedef.h"
#include "data.h"
#include "stream-tools.h"

/// Because of the god awful mess that are vcf header lines.
/** This is likely to become some form of container to handel moving data into and out of rows of map file.
 */

class Vcf_header : public virtual Data_interface {
private:
public:

};

#endif
class Vcf_data : public virtual Data_interface {
private:
	using Indexed_data::abs_pos_;
public:
	std::vector <std::string> ID;
	Base ref;
	std::vector <std::string> ALT;
	float_t QUAL;
	bool FILTER;
	std::vector <std::string> INFO;

        AA ancestral allele
        AC allele count in genotypes, for each ALT allele, in the same order as listed
        AF allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
        AN total number of alleles in called genotypes
        BQ RMS base quality at this position
        CIGAR cigar string describing how to align an alternate allele to the reference allele
        DB dbSNP membership
        DP combined depth across samples, e.g. DP=154
        END end position of the variant described in this record (esp. for CNVs)
        H2 membership in hapmap2
        MQ RMS mapping quality, e.g. MQ=52
        MQ0 Number of MAPQ == 0 reads covering this record
        NS Number of samples with data
        SB strand bias at this position
        SOMATIC indicates that the record is a somatic mutation, for cancer genomics
        VALIDATED validated by follow-up experiment

};

#endif
