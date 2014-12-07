/* 

Program Estimator.cpp:

	1) input from a list of site-specific quartets from a pooled population sample;

	2) identify the major and minor alleles, obtain the maximum-likelihood estimates of allele frequencies, and significance levels.

The designated major nucleotide is simply the one with the highest rank, and the minor nucleotide the one with the second highest rank.
	If the top three ranks are all equal, the site is treated as unresolvable, with both major and minor nucelotides designated in the output by a *.
	If the second and third ranks are equal but lower than the major-nucleotide count, the site is treated as monomorphic, with the minor nucleotide designated by a *.
	
Input File: six columns; first two are arbitrary identifiers of the site (e.g., chromosome and position); final four are integer numbers of reads of A, C, G, T.
	Columns are tab delimited.
	Default name is "datain.txt".

Output File: two columns of site identifiers; major allele; minor allele; major-allele frequency; minor-allele frequency; error rate; likeilhood-ratio test of polymorphism. 
	Columns are comma delimited.
	Default name is "dataout.txt".

***** Note that significance at the 0.05, 0.01, 0.001 levels requires that the likelihood-ratio test statistic exceed 3.841, 6.635, and 10.827, respectively. 

***** The 95% support interval can be obtained by determining the changes in the estimate of p in both directions required to reduce the log likelihood by 3.841
      although this is not currently implemented. 
*/

#include "Estimator.hpp"

void estimator(const char *infile, const char *outfile)
{

/* Point to the input and output files. */

FILE *instream=NULL;
FILE *outstream=NULL;

int n[5];						/* read quartet for numbers of reads; 1 = A; 2 = C; 3 = G; 4 = T. */

int coverage;						/* total number of reads at the site */

int nmax;						/* number of reads for the most abundant type (putative major allele) */
int nmin;						/* number of reads for the second most abundant type (putative minor allele) */
int ner1;						/* number of reads for the third highest ranked nucleotide */

int major, minor;					/* identities of major and minor alleles */
int third;						/* identity of the third-place nucleotide */

double pmaj;						/* fraction of putative major reads among the total of major and minor */

double eml, pml;					/* maximum-likelihood estimates of the error rate and major-allele frequency */

char nuc[6];						/* converts numerical to alphabetical indicator of nucleotide */

char c;							/* used to peek at the next character in the file */

double lphimaj, lphimin, lphie;				/* terms for the likelihood expression */
double emltheta;

double llhoodp, llhoodm;				/* log likelihoods under the assumption of polymorphism and monomorphism */

int nerrors, nerrorsm;					/* errors inferred in the full- and reduced-model likelihood analyses */

double errorm;						/* estimated error rate in the reduced model */

double llstat;						/* difference in the log likelihoods */

char id1[30];						/* identifiers for the site in the first two columns of the input file */
char id2[30];


nuc[1] = 'A';
nuc[2] = 'C';
nuc[3] = 'G';
nuc[4] = 'T';
nuc[5] = '*';

/* Open the output and input files. */

outstream = fopen(outfile, "w");
instream = fopen(infile, "r");
if(outstream==NULL) {std::cerr << "failed to open file " << outfile << " for writing\n"; exit(0); }
if(instream==NULL) {std::cerr << "failed to open file " << infile << " for reading\n"; exit(0); }

/* Start inputting and analyzing the data line by line. */

while (! feof (instream) ){

c = fgetc (instream);
fseek(instream, -1, SEEK_CUR);

/* check for start of new scaffold*/

if (c!='>'){if (fscanf(instream,"%s\t%i\t%i\t%i\t%i", id2, &n[1], &n[2], &n[3], &n[4])==EOF) break;}
else {
	if (fscanf(instream, "%s", id1)==EOF) break;
	continue;
	};


coverage = n[1] + n[2] + n[3] + n[4];

/* Identify the counts for the putative major and minor alleles.  */

if (n[1] >= n[2]) {
	nmax = n[1]; major = 1;
	nmin = n[2]; minor = 2; } 
else {
	nmax = n[2]; major = 2;
	nmin = n[1]; minor = 1; }

if (n[3] > nmax) {
	ner1 = nmin; third = minor;
	nmin = nmax; minor = major;
	nmax = n[3]; major = 3; } 
else if (n[3] > nmin) {
	ner1 = nmin; third = minor;
	nmin = n[3]; minor = 3; }
else {ner1 = n[3]; third = 3;}

if (n[4] > nmax) {
	ner1 = nmin; third = minor;
	nmin = nmax; minor = major;
	nmax = n[4]; major = 4;
} else if (n[4] > nmin) {
	ner1 = nmin; third = minor;
	nmin = n[4]; minor = 4; }
else if (n[4] > ner1) {
	ner1 = n[4]; third = 4; }
	
/* Calculate the ML estimates of the major / minor allele frequencies and the error rate. */

pmaj = ((double) nmax) / ( ((double) nmax) + ((double) nmin) );

eml = 1.5 * ( ((double) coverage) - ((double) nmax) - ((double) nmin) ) / ((double) coverage);
	
pml = (pmaj * (1.0 - (2.0 * eml / 3.0))) - (eml / 3.0);
pml = pml / (1.0 - (4.0 * eml / 3.0));

/* Calculate the log likelihoods under: 1) the full model; and 2) the reduced model assuming no polymorphism (major-allele frequency = 1.0). */

/* First, get the log likelihood under the full model, which assumes a polymorphism. */
	
emltheta = 1.0 - (4.0 * eml / 3.0);
lphimaj = (pml * emltheta) + (eml / 3.0);
lphimin = 1.0 - eml - (pml * emltheta);
lphie =  eml / 3.0;
		
nerrors = coverage - nmax - nmin;

if ( (lphie == 0.0) && (lphimin > 0.0) ) {
	llhoodp = (((double) nmax) * log(lphimaj)) + (((double) nmin) * log(lphimin));}
else if ( (lphie == 0.0) && (lphimin == 0.0) ) {
	llhoodp = 0.0;}
else {
	llhoodp = (((double) nmax) * log(lphimaj)) + (((double) nmin) * log(lphimin)) + (((double) nerrors) * log(lphie));}
	

/* Second, get the log likelihood under the assumption of just one true allele. */

nerrorsm = coverage - nmax;
		
errorm = ((double) nerrorsm) / ((double) coverage);
		
if ( errorm == 0.0 ) {
	llhoodm = 0.0;}
else {
	llhoodm = ( ((double) nmax) * log(1.0 - errorm) ) + (((double) nerrorsm) * log(errorm/3.0) );}


/* Likelihood ratio test statistic; asymptotically chi-square distributed with one degree of freedom. */

llstat = 2.0 * (llhoodp - llhoodm);



/* If the top three ranks are all equal, flag as unresolvable by denoting the alleles with asterisks. */

if  ( (nmax == nmin) && (nmin == ner1) ) {
	pml = 0.0;
	eml = 0.0;
	llstat = 0.0; 
	major = 5;
	minor = 5;}


/* If the second and third ranks are equal but less than the top, treat as fixed for the putative major allele. */

if  ( (nmax > nmin) && (nmin == ner1) ) {
	pml = 1.0;
	llstat = 0.0; 
	minor = 5;}



printf("%s\t%s\t%c\t%c\t%6.5f\t%6.5f\t%6.5f\t%5d\t%10.5f\n", id1 , id2 , nuc[major] , nuc[minor] , pml , (1.0 - pml) , eml , coverage, llstat); 
fprintf(outstream, "%s\t%s\t%c\t%c\t%6.5f\t%6.5f\t%6.5f\t%5d\t%10.5f\n", id1 , id2 , nuc[major] , nuc[minor] , pml , (1.0 - pml) , eml , coverage, llstat); 


}					/* Ends the computations when the end of file is reached. */

printf("\n");

exit(0);

}


