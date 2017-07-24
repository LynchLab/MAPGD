#include "vcf-file.h"

#ifndef NOHTS

void call_minor(const  Genotype &gtl, int32_t *calls)
{
	if (gtl.MM>gtl.Mm)
	{
		if(gtl.MM>gtl.mm)
		{
			calls[0]=bcf_gt_unphased(1);
			calls[1]=bcf_gt_unphased(1);
		} else {
			calls[0]=bcf_gt_unphased(0);
			calls[1]=bcf_gt_unphased(0);
		}
	} else if(gtl.Mm>gtl.mm) {
		calls[0]=bcf_gt_unphased(0);
		calls[1]=bcf_gt_unphased(1);
	} else {
		calls[0]=bcf_gt_unphased(0);
		calls[1]=bcf_gt_unphased(0);
	}
}

void call_major(const  Genotype &gtl, int32_t *calls)
{
	if (gtl.MM>gtl.Mm)
	{
		if(gtl.MM>gtl.mm)
		{
			calls[0]=bcf_gt_unphased(0);
			calls[1]=bcf_gt_unphased(0);
		} else {
			calls[0]=bcf_gt_unphased(1);
			calls[1]=bcf_gt_unphased(1);
		}
	} else if(gtl.Mm>gtl.mm) {
		calls[0]=bcf_gt_unphased(0);
		calls[1]=bcf_gt_unphased(1);
	} else {
		calls[0]=bcf_gt_unphased(1);
		calls[1]=bcf_gt_unphased(1);
	}
}

Vcf_data::Vcf_data ()
{
	record_=NULL;
}

void 
Vcf_file::read (Vcf_data &vcf)
{
	if (open_) 
	{
		std::cerr << __FILE__ << ":" << __LINE__ <<  "trying to read.\n";
		if (bcf_read(file_, vcf.header_, vcf.record_)!=0) 
		{
			std::cerr << __FILE__ << ":" << __LINE__ <<  "attempt to read from file failed.\n";
			exit(0);
		}
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ <<  "attempt to read from unopened file.\n";
		exit(0);
	}
}

void
Vcf_data::put(const File_index &index, const Allele &allele, const Population &pop)
{
	if(allele.get_abs_pos()!=pop.get_abs_pos() )
	{
		std::cerr << __FILE__ << ":" << __LINE__ <<  "allele and locus have different positions (private abs_pos_).\n";
	}

	char alleles[4]={0};
	size_t size=record_->n_sample;
	id1_t abs_pos=pop.get_abs_pos();
	float *gp=new float[size*3], freq=allele.freq;
	float *vit=gp;
	int32_t *dp=new int32_t[size];
	int32_t *dp_it=dp;
	int32_t *gt=new int[size*2];
	int32_t *gt_it=gt;
	record_->rid = index.get_id0(abs_pos);
	record_->pos = index.get_id1(abs_pos)-1;
	if (allele.major==allele.ref)
	{
		sprintf(alleles,"%c,%c", Base::btoc(allele.ref), Base::btoc(allele.minor) );

        	bcf_update_info_float(header_, record_, "AF", &freq, 1);
       // 	bcf_update_info_float(header_, record_, "", &f, 1);

		for (std::vector<Genotype>::const_iterator it=pop.likelihoods.cbegin(); it<pop.likelihoods.cend(); ++it)
       		{
			*(vit+0)=(float)log(it->MM);
			*(vit+1)=(float)log(it->Mm);
			*(vit+2)=(float)log(it->mm);
			*dp_it=(int32_t)(it->N);
			call_major(*it, gt_it);
			vit+=3;
			++dp_it;
			gt_it+=2;
		}
	} else {
		freq=1.-freq;
		sprintf(alleles,"%c,%c", Base::btoc(allele.ref), Base::btoc(allele.major) );

        	bcf_update_info_float(header_, record_, "AF", &freq, 1);
     //   	bcf_update_info_float(header_, record_, "", &f, 1);

		for (std::vector<Genotype>::const_iterator it=pop.likelihoods.cbegin(); it<pop.likelihoods.cend(); ++it)
		{
			*(vit+0)=(float)log(it->mm);
			*(vit+1)=(float)log(it->Mm);
			*(vit+2)=(float)log(it->MM);
			*dp_it=(int32_t)(it->N);
			call_minor(*it, gt_it);
			vit+=3;
			++dp_it;
			gt_it+=2;
		}
	}
	bcf_update_alleles_str(header_, record_, alleles);
        bcf_update_format_float(header_, record_, "GP", gp, size*3);
	bcf_update_genotypes(header_, record_, gt, size*2);
        bcf_update_format_int32(header_, record_, "DP", dp, size);
	free(gp);
	free(gt);
	free(dp);
}

void
Vcf_data::put(const Data *data, ...)
{
	va_list args;
	va_start(args, data);
	File_index *idx = dynamic_cast<File_index *>(va_arg(args, Data *) );
	Population *pop = dynamic_cast<Population *>(va_arg(args, Data *) );

//	put(*idx, *pop);
}

id1_t
Vcf_data::get(const File_index &index, Population &pop) const
{
	char alleles[4];
	size_t size=record_->n_sample;


	float *genotype_qualities=new float[sample_size_*3];
	float *vit=genotype_qualities;
	int *sequence_depth=new int[sample_size_];
	int *dit=sequence_depth;	

	int genotype_qualities_size=sample_size_*3;
	int sequence_depth_size=sample_size_;

	bcf_get_format_float(header_, record_, "GQ", &genotype_qualities, &genotype_qualities_size);
	bcf_get_format_int32(header_, record_, "DP", &sequence_depth, &sequence_depth_size);

	std::vector<Genotype>::iterator end=pop.likelihoods.end();
	for (std::vector<Genotype>::iterator it=pop.likelihoods.begin(); it<end; ++it)
	{
			it->mm=(float)(exp(vit[0]) );
			it->Mm=(float)(exp(vit[1]) );
			it->MM=(float)(exp(vit[2]) );
			it->N=*(++dit);
			vit+=3;
	}

	int allele_count=1;
	bcf_get_info_float(header_, record_, "AF", &pop.m, &allele_count);
	//bcf_get_info_float(header_, record_, "", &pop.f, 1);

	if (pop.m>0.5)
	{
		pop.major=Base::ctob(record_->d.allele[0][0]);
		pop.minor=Base::ctob(record_->d.allele[1][0]);
	}
	else if (pop.m < 0.5)
	{
		pop.major=Base::ctob(record_->d.allele[0][0]);
		pop.minor=Base::ctob(record_->d.allele[1][0]);
	}
	else 
	{
		pop.major=Base::ctob(record_->d.allele[0][0]);
		pop.minor=Base::ctob(record_->d.allele[1][0]);
	}
	pop.set_abs_pos(index.get_abs_pos(record_->rid, record_->pos+1) );
	return 0;
}


void 
Vcf_data::get(Data *data, ...) const
{
	if (0)
	{
	//	get_text("GOOD WORK %s", "BOO!");
	} else {
	//	get_text("BAD WORK %s", "BOO!");
	}
}

//Takes a vcf file and prints a bcf to stream.
void 
Vcf_file::write(const Vcf_data &vcf)
{
	bcf_write(file_, vcf.header_, vcf.record_);
}

void
Vcf_file::write_header(const Vcf_data &vcf)
{
	bcf_hdr_write(file_, vcf.header_);
}

void
Vcf_file::open(const std::ios_base::openmode &mode)
{
//	if (mode & READ) open(std::cin, mode); 
//	else if (mode & WRITE) open(std::cout, mode);
//	else {
		std::cerr << __FILE__ << ":" << __LINE__ << ": cannot open file in std::ios_base::openmode &mode=" << mode << std::endl; 
//	}
}


/*
void 
Vcf_file::close(void)
{
}*/

void
Vcf_data::set_header(const File_index &index, const std::vector <std::string> &sample_names)
{

/* === Dictionary ===

   The header keeps three dictionaries. The first keeps IDs in the
   "FILTER/INFO/FORMAT" lines, 
   The second keeps the sequence names and lengths in the "contig" lines 
   popultation.get_index() and the last keeps the sample names 
   get_header().get_sample_names();

   bcf_hdr_t::dict[]

   Is the actual hash table, which is opaque to the end users. In the hash
   table, the key is the ID or sample name as a C string and the value is a
   bcf_idinfo_t struct. bcf_hdr_t::id[] points to key-value pairs in the hash
   table in the order that they appear in the VCF header. bcf_hdr_t::n[] is the
   size of the hash table or, equivalently, the length of the id[] arrays.
*/
	
	header_=bcf_hdr_init("w");

	bcf_hdr_append(header_, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">");
	bcf_hdr_append(header_, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">");
	bcf_hdr_append(header_, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	bcf_hdr_append(header_, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Probabilities\">");
	bcf_hdr_append(header_, "##FORMAT=<ID=DP,Number=1,Type=Float,Description=\"Read Depth\">");
//	if (args->record_cmd_line) 
//	bcf_hdr_append_version(args->header, args->argc, args->argv, "mapgd convert");
	for (std::vector <std::string>::const_iterator sample=sample_names.cbegin(); sample!=sample_names.cend(); ++sample )
	{
		if (bcf_hdr_add_sample(header_, sample->c_str() ) );
	}
	std::vector<std::string> scaffold_names=index.get_names();
	std::vector<id1_t> scaffold_sizes=index.get_sizes();

	for(size_t x=0; x<scaffold_names.size(); ++x)
	{
		bcf_hdr_printf(header_, "##contig=<ID=%s,length=%d>", scaffold_names[x].c_str(), scaffold_sizes[x]);   // MAX_CSI_COOR	
	}
	//if(!record_) 
	record_=bcf_init();
	bcf_hdr_sync(header_);
	record_->n_sample=bcf_hdr_nsamples(header_);
//
}

void
Vcf_file::open(const char *file_name, const std::ios_base::openmode &mode)
{
	if (mode & READ)
	{
		std::cerr << "opening\n";
		file_=bcf_open(file_name, "r");
		if ( file_) open_=true;
		else {
			std::cerr << __FILE__ << __LINE__ << "Whoops! bad bad bad\n";
			exit(EXIT_FAILURE);
		}
	}
	else if (mode & WRITE)
	{
		file_=vcf_open(file_name, "w");
		if ( file_) open_=true;
	}
	if (open_)
	{
		if (mode & READ)
		{
			// Read header
			// ALOCATE HEADER and BCF RECORD
			// bcf_init
		}
	} else {
		std::cerr << __FILE__ << ":" << __LINE__ << " failure to open file.\n";
	}
		
}
void 
Vcf_file::close(void)
{
	if (open_)
	{
		vcf_close(file_);
		open_=false;
	}
		
}

File_index
Vcf_data::get_index (void) const
{
	File_index index;
	int id_size=header_->n[BCF_DT_CTG];
        for (int i = 0; i < id_size; i++) {
		index.add_id(header_->id[BCF_DT_CTG][i].key, header_->id[BCF_DT_CTG][i].val->info[0] ); 
        }
	return index;
}
	
std::vector<std::string> 
Vcf_data::get_sample_names (void) const 
{
	std::vector<std::string> names;
	int sample_size=header_->n[BCF_DT_SAMPLE];
	for (int i = 0; i < sample_size; i++) {
		names.push_back(header_->id[BCF_DT_SAMPLE][i].key );
	}       
	return names;
}

Vcf_data
Vcf_file::read_header(void)
{
	std::cerr << "Reading header dumbass\n";
	Vcf_data vcf;
	vcf.header_=bcf_hdr_read(file_);
	vcf.record_=bcf_init();
	vcf.record_->n_sample=bcf_hdr_nsamples(vcf.header_);
	vcf.sample_size_=vcf.record_->n_sample;
//	vcf.max_unpack=?
	table_open_=true;
	return vcf;
}
#endif
