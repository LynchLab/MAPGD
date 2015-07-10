/* 
Sub-routine: PooledComparison.cpp ver_1: 

1) use ML to identify the major and minor alleles across populations, estimate the major-allele frequencies, in retrospect.

2) test for the significance of allele-frequency differences.
	
It is assumed that data are for multiple diploid individuals from populations in Hardy-Weinberg equilibrium, which is 
	equivalent to 2*nsample haploid sets of chromosomes being sampled.

Sequencing errors are assumed to be random and equal in all directions.

Note that in the analysis, the two most abundantly read nucleotides are interpreted as the major and minor alleles in the pooled population.
	When the sample size is small, and/or the true frequencies are near 0.5, this means that the actual major and minor
	alleles may be misidentified. 

The error rate is estimated from the pooled sample.

For more information please read 

********** Within the program, just below the list of terms, the set of minor-allele frequencies that one wishes to compute for one population can be set. **********
********** The allele frequency of the second population is entered at the top as a #define statement. **********

*/

/* ********************************************************************************************************** */

#define niters		50000		/* number of iterations of the sampling / estimation procedure */
#define EMLMIN		0.000001	/* minimum error estimate				       */

#include	<stdio.h>
#include 	<math.h>
#include 	<sys/types.h>
#include	<stdlib.h>
#include	<time.h>
#include 	"ReadFile.hpp"

/* ********************************************************************************************************** */

/* Point to the output file. */

void compare(const char *infile1, const char *infile2, const char *outfile)
{

FILE *output;

/* ************************************************************************************************************ */

int major, minor;					/* indicators for observed major and minor alleles */

double minfreq[100];					/* minor-allele frequencies at which the computations will be done */
int nminor;						/* counter for loop of minor-allele frequencies */

double majfreq;						/* major-allele frequency */

double hetero;						/* heterozygosity under the assumption of Hardy-Weinberg equilibrium */
double hommaj;						/* homozygote frequencies under Hardy-Weinberg */

double hetcut;						/* cutoff for heterozygotes in random sampling from population */

double popsamp;						/* population sample size */

double majfreqB, heteroB, hommajB, hetcutB;

double popin;						/* individual randomly drawn from population */
double seqin;						/* read randomly drawn from the sample */

double etheta, mincut1, mincut2;	/* constants used to compute the draws of reads from the population sample */
double mlterm;

int n1read[5];						/* numbers of the four nucleotide reads for the actual major and minor alleles and the two error types */
int n2read[5];						/* for populations 1, 2 , and pooled */
int n3read[5];	

double freqmaj;						/* frequency of the major allele in the sample */

double phimaj;						/* expected frequency of draws of the major allele after accounting for errors */

int nmaxP, nminP;					/* counts for the putative major and minor alleles in the pooled sample */

double pmaj;					/* fraction of putative major reads among the total of major and minor */
double pmajA, pmajB;

double eml, pml;				/* maximum-likelihood estimates of the error rate and major-allele frequency */
double pmlA, pmlB;

double totests;						/* sums to use in estimating mean and standard deviation of the estimates over all runs */

double emltheta;

double coef;					/* binomial coefficents used in the test statistic */

double lnterm;
double lnpml, ln1pml;
//WTF?!?!?!?
double exphimaj[20001];				/* terms used in the likelihood test statistic */
double exphimin[20001];

double lhoodP1, lhoodS1;			/* likelihoods of the data in population 1 under the assumptions of pooled or distinct allele frequencies */
double lhoodP2, lhoodS2;			/* likelihoods of the data in population 2 under the assumptions of pooled or distinct allele frequencies */

double llhoodP, llhoodS;			/* log likelihoods under the models of homogeneity and heterogeneity of allele frequencies */ 

double var1, sdsamp1;				/* theoretical variance and SD of allele frequency */
double var2, sdsamp2;	

double minsearch, minsearch1, minsearch2;		/* cutoff point indicators for the likelihood function, beyond which there should be no significant contribution */
double maxsearch, maxsearch1, maxsearch2;

int minig, maxig, ig;					/* limits to allele frequency count used in the likelihood analysis, to streamline the computational time */

int nsample;
int n2samp;
int coverage;
int minfreq2;

double llstat;						/* difference in the log likelihoods */
double sig1, sig2, sig3, sig4;			/* counters for significance tests */

double lnfact[2010];

profile pop[2];
/* Open the output file. */

pop[0].open(infile1, "r");
pop[1].open(infile2, "r");
output=fopen(outfile, "w");

/* ********** Set the minor-allele frequencies for the first population at which the computations will be done. ************* */
/*
Set the initial counters for the statistics of the simulations equal to zero. */

totests = 0.0;
sig1 = 0.0;
sig2 = 0.0;
sig3 = 0.0;
sig4 = 0.0;

/* Calculate the starting and ending points for the allele counts used in the likelihood function. */

/*
else { maxig = int(maxsearch); }
*/

minig = 0;
maxig = 2 * nsample;

/* 1) GET THE READS FOR THE FIRST POPULATION. */
/* 2) GENERATE THE READS FOR THE SECOND POPULATION. */

	sync(2, pop)

/* 3) GENERATE THE READ COUNTS POOLED OVER BOTH SAMPLES */

	n3read[1] = n1read[1] + n2read[1];
	n3read[2] = n1read[2] + n2read[2];
	n3read[3] = n1read[3] + n2read[3];
	n3read[4] = n1read[4] + n2read[4];

/* 4) IDENTIFY THE OVERALL MAJOR AND MINOR ALLELES. */ 

	/* Identify the counts for the Putative major and minor alleles	for the pooled pair of populations. */ 

	if (n3read[1] > n3read[2]) {
		nmaxP = n3read[1]; major = 1;
		nminP = n3read[2]; minor = 2;} 
	else {
		nmaxP = n3read[2]; major = 2;
		nminP = n3read[1]; minor = 1;}
	
	if (n3read[3] > nmaxP) {
		nminP = n3read[major]; minor = major;
		nmaxP = n3read[3]; major = 3;} 
	else if (n3read[3] > nminP) {
		nminP = n3read[3]; minor = 3;}
	
	if (n3read[4] > nmaxP) {
		nminP = n3read[major]; minor = major;
		nmaxP = n3read[4]; major = 4;} 
	else if (n3read[4] > nminP) {
		nminP = n3read[4]; minor = 4;}


/* 5) CALCULATE THE LIKELIHOOD OF THE JOINTLY POOLED POPULATIONS. */
	
	/* Calculate the ML estimates of the major pooled allele frequency and the error rate. */

	pmaj = ((double) n3read[major]) / ( ((double) n3read[major]) + ((double) n3read[minor]) );
	eml = 1.5 * ( (2.0 * (double) coverage) - ((double) n3read[major]) - ((double) n3read[minor]) ) / (2.0 * (double) coverage);
		if (eml == 0.0) { eml = EMLMIN;}
	pml = (pmaj * (1.0 - (2.0 * eml / 3.0))) - (eml / 3.0);
	pml = pml / (1.0 - (4.0 * eml / 3.0));
	emltheta = 1.0 - (4.0 * eml / 3.0);
	
	
	/* Calculate the likelihoods under the reduced model assuming no variation between populations. */
			
	/* Calculate the component for the first population. */

		lhoodP1 = 0.0;
		
		for (ig = minig; ig <= maxig; ++ig) {
			exphimaj[ig] = log(((((double) ig) / popsamp) * emltheta) + (eml / 3.0));
			exphimin[ig] = log(1.0 - eml - ((((double) ig) / popsamp) * emltheta)); }

		if ( (pml *(1.0-pml)) > 0.0) {
		
			lnpml = log(pml);
			ln1pml = log(1.0 - pml);

			for (ig = minig; ig <= maxig; ++ig) {
				lnterm = lnfact[2*nsample] - lnfact[ig] - lnfact[(2*nsample)-ig] + (((double) ig) * lnpml) + (((double) (n2samp-ig)) * ln1pml)
					+ (((double) n1read[major]) * exphimaj[ig]) +  (((double) n1read[minor]) * exphimin[ig]); 
				lhoodP1 = lhoodP1 + exp(lnterm); }
		}
		
		/* Calculate the component for the second population. */

		lhoodP2 = 0.0;
		if ( (pml *(1.0-pml)) > 0.0) {
			for (ig = minig; ig <= maxig; ++ig) {
				lnterm = lnfact[2*nsample] - lnfact[ig] - lnfact[(2*nsample)-ig] + (((double) ig) * lnpml) + (((double) (n2samp-ig)) * ln1pml)
					+ (((double) n2read[major]) * exphimaj[ig]) +  (((double) n2read[minor]) * exphimin[ig]); 
				lhoodP2 = lhoodP2 + exp(lnterm); }
		}
		
		if ( (pml *(1.0-pml)) <= 0.0) {
			llhoodP = 0.0;}
		else { llhoodP = log(lhoodP1) + log(lhoodP2); }

		/* 6) CALCULATE THE LIKELIHOOD UNDER THE ASSUMPTION OF POPULATION SUBDIVISION. */ 

		pmajA = ((double) n1read[major]) / ( ((double) n1read[major]) + ((double) n1read[minor]) );
		pmajB = ((double) n2read[major]) / ( ((double) n2read[major]) + ((double) n2read[minor]) );
	
		pmlA = (pmajA * (1.0 - (2.0 * eml / 3.0))) - (eml / 3.0);
		pmlA = pmlA / (1.0 - (4.0 * eml / 3.0));
		pmlB = (pmajB * (1.0 - (2.0 * eml / 3.0))) - (eml / 3.0);
		pmlB = pmlB / (1.0 - (4.0 * eml / 3.0)); 

		if (pmlA >= 1.0) { pmlA = 0.999999999; }
		if (pmlB >= 1.0) { pmlB = 0.999999999; }

		/* Calculate the component for the first population. */

		lhoodS1 = 0.0;
		lnpml = log(pmlA);
		ln1pml = log(1.0 - pmlA);
	
		if (  ((pmlA * (1.0 - pmlA) ) != 0.0)  ) { 
			for (ig = minig; ig <= maxig; ++ig) {
				lnterm = lnfact[2*nsample] - lnfact[ig] - lnfact[(2*nsample)-ig] + (((double) ig) * lnpml) + (((double) (n2samp-ig)) * ln1pml)
					+ ( ((double) n1read[major]) * exphimaj[ig]) +  (((double) n1read[minor]) * exphimin[ig]); 
				lhoodS1 = lhoodS1 + exp(lnterm); } } 
		
		/* Calculate the component for the second population. */

		lhoodS2 = 0.0;
		lnpml = log(pmlB);
		ln1pml = log(1.0 - pmlB);
		
		if ( ((pmlB * (1.0 - pmlB) ) != 0.0) ) { 
			for (ig = minig; ig <= maxig; ++ig) {
				lnterm = lnfact[2*nsample] - lnfact[ig] - lnfact[(2*nsample)-ig] + (((double) ig) * lnpml) + (((double) (n2samp-ig)) * ln1pml)
					+ ( ((double) n2read[major]) * exphimaj[ig]) +  (((double) n2read[minor]) * exphimin[ig]); 
				lhoodS2 = lhoodS2 + exp(lnterm); } }

		if ( ((lhoodS1 * lhoodS2) == 0.0) ) {
			llhoodS = 0.0;}

		else { llhoodS = log(lhoodS1) + log(lhoodS2);} 

		
	/* Likelihood ratio test statistic; asymptotically chi-square distributed with one degree of freedom. */

	llstat = 2.0 * (llhoodS - llhoodP);

	/* Significance at the 0.05, 0.01, 0.001 levels requires the statistic, with 1 degrees of freedom, to exceed 3.841, 6.635, and 10.827, respectively. */

	if (llstat > 3.841) {
		sig1 = sig1 + 1.0;}
	if (llstat > 6.635) {
		sig2 = sig2 + 1.0;}
	if (llstat > 10.828) {
		sig3 = sig3 + 1.0;}
	if (llstat > 15.14) {
		sig4 = sig4 + 1.0;}

	totests = totests + 1.0;

}	/* End the loop for population samples for this specific minor-allele frequency in the first population. */


/* Calculate the means and variances of estimates over all individual runs for this minor-allele frequency. */
printf("Sign. at 0.05, 0.01, 0.001, 0.0001 levels = %5.4f , %5.4f , %5.4f , %5.4f\n\n", (sig1/totests) , (sig2/totests) , (sig3/totests), (sig4/totests));

fprintf(output, "%8.7f, %d, %d, %d, %8.7f, %8.7f, %8.7f, %8.7f\n",
		eml, coverage, nsample, niters, (sig1/totests), (sig2/totests), (sig3/totests), (sig4/totests) );

exit(0);

}
