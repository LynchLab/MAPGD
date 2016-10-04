#include "vcf-file.h"

#ifndef NOHTS

Vcf_data::Vcf_data ()
{
}

void Vcf_data::read (Allele &A) 
{
}

void Vcf_data::write (Allele &A) const
{
/*	char alleles[2][1];
	alleles[0][0]='T';
	alleles[0][1]='A';
	bcf_update_alleles(header, record, alleles, 2);
*/	//bcf_update_info(bcf_hdr_t *hdr, bcf1_t *line, const char *key, const void *values, int n, int type)
//	bcf_write(std::cout, header, record);
}

void Vcf_file::read (Allele &A) 
{
	if (open_) 
	{
		std::string line;
		std::getline(*in_, line);
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ << " attempt to write to unoppend file.\n";
		exit(0);
	}
}

void Vcf_file::write (const Allele &A) const
{
	if (open_)
	{
		std::cerr << "is open.\n";
		std::string line="doot!\n";
		std::cerr << line;
		//vcf_data.write(*out_, A);
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ << " attempt to write to unoppend file.\n";
		exit(0);
	}
}

#if 0
Vcf_file::open(const std::ios_base::openmode &mode)
{
}

Vcf_file::set_header(?, Indexed_data <Population> &population)
{

/* === Dictionary ===

   The header keeps three dictonaries. The first keeps IDs in the
   "FILTER/INFO/FORMAT" lines, 

   the second keeps the sequence names and lengths in the "contig" lines 

   popultation.get_index()

   and the last keeps the sample names 

   get_header().get_sample_names();

   bcf_hdr_t::dict[]

   is the actual hash table, which is opaque to the end users. In the hash
   table, the key is the ID or sample name as a C string and the value is a
   bcf_idinfo_t struct. bcf_hdr_t::id[] points to key-value pairs in the hash
   table in the order that they appear in the VCF header. bcf_hdr_t::n[] is the
   size of the hash table or, equivalently, the length of the id[] arrays.
*/
	if (open_)
	{
		header=bcf_hdr_init();
		while (std::vector <std::string>::const_iterator sample=samples.cbegin(); ++sample; sample!=samples.cend() )
		{
			 if (bcf_hdr_add_sample(header, sample->c_str() ) );
		}
		// const char **bcf_hdr_seqnames(header, int *nseqs)
		
	// ALOCATE HEADER and BCF RECORD
	// bcf_init
	}
}

Vcf_file::open(const char *file_name, const std::ios_base::openmode &mode)
{
	open_no_extention(file_name, mode);
	if (open_)
	{
		if (mode && READ)
		{
			// Read header
			// ALOCATE HEADER and BCF RECORD
			// bcf_init
		}
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ << " failure to open file.\n";
	}
		
}
#endif

#endif
