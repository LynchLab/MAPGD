#include "simple_types.h"

GENPROB::GENPROB(void)
{
	keynum_=GENPROB;
	keyname_=16;
	keydesc_=genotypic probability;
}
FREQ::FREQ(void)
{
	keynum_=FREQ;
	keyname_=17;
	keydesc_=allele frequency;
}
ERROR::ERROR(void)
{
	keynum_=ERROR;
	keyname_=18;
	keydesc_=error rate;
}
LOGLIKE::LOGLIKE(void)
{
	keynum_=LOGLIKE;
	keyname_=19;
	keydesc_=log likelihood;
}
SMPNAME::SMPNAME(void)
{
	keynum_=SMPNAME;
	keyname_=20;
	keydesc_=sample name;
}
ROWID::ROWID(void)
{
	keynum_=ROWID;
	keyname_=21;
	keydesc_=a number representing how many rows occur before this row in the file;
}
COV::COV(void)
{
	keynum_=COV;
	keyname_=22;
	keydesc_=read depth;
}
SMPSIZE::SMPSIZE(void)
{
	keynum_=SMPSIZE;
	keyname_=23;
	keydesc_=the number of samples/individuals within the data;
}
SOMATIC::SOMATIC(void)
{
	keynum_=SOMATIC;
	keyname_=24;
	keydesc_=a flag to indicate that the record is a somatic mutaiton;
}
VALID::VALID(void)
{
	keynum_=VALID;
	keyname_=25;
	keydesc_=a flag to indicate that a record has been validated?;
}
STRBIAS::STRBIAS(void)
{
	keynum_=STRBIAS;
	keyname_=26;
	keydesc_=strand bais at this position;
}
SB::SB(void)
{
	keynum_=SB;
	keyname_=26;
	keydesc_=strand bais at this position;
}
MAPQ0::MAPQ0(void)
{
	keynum_=MAPQ0;
	keyname_=27;
	keydesc_=number of MAPQ==0 reads covering this position;
}
HAPMAP2::HAPMAP2(void)
{
	keynum_=HAPMAP2;
	keyname_=28;
	keydesc_=membership in hapmap2;
}
END::END(void)
{
	keynum_=END;
	keyname_=29;
	keydesc_=ending row of a variant described by this record;
}
REF::REF(void)
{
	keynum_=REF;
	keyname_=30;
	keydesc_=reference allele;
}
ANCTYPE::ANCTYPE(void)
{
	keynum_=ANCTYPE;
	keyname_=31;
	keydesc_=ancestoral allele;
}
NCUT::NCUT(void)
{
	keynum_=NCUT;
	keyname_=32;
	keydesc_=the number of samples excluded because of filters;
}
GOF::GOF(void)
{
	keynum_=GOF;
	keyname_=33;
	keydesc_=goodness of fit value;
}
DEFAULT::DEFAULT(void)
{
	keynum_=DEFAULT;
	keyname_=34;
	keydesc_=a flag to indicate that a record represents a default value;
}
POLYLL::POLYLL(void)
{
	keynum_=POLYLL;
	keyname_=35;
	keydesc_=the log of the ratio polymorphic fit / monomorphic fit (chi-sqr, 2 df);
}
HWELL::HWELL(void)
{
	keynum_=HWELL;
	keyname_=36;
	keydesc_=the log of the ratio polymorphic fit / monomorphic fit (chi-sqr, 2 df);
}
MAXLL::MAXLL(void)
{
	keynum_=MAXLL;
	keyname_=37;
	keydesc_=the maximum likelihood fit of the data;
}
EFCHROM::EFCHROM(void)
{
	keynum_=EFCHROM;
	keyname_=38;
	keydesc_=the effective number of chromosomes at a locus;
}
