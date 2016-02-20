#MAPGD version 0.4.2

#Contents
#Introduction
#Basic Design
#Tutorials
#Introduction

MAPGD is a series of related programs that estimate allele frequency, heterozygosity, Hardy-Weinberg disequilibrium and identity by descent (IBD) coefficients from population genomic data using statistically rigorous maximum-likelihood approach.
Basic Design:
MAPGD:

#Style Guidelines

mapgd is written to conform to the GNU style guidelines, at least to the extent that I have had time to read and implement the guidelines.

* Class names should be whole English words starting with an initial capital (to distinguish them from the core classes and types of c++): e.g My\_class\_name.
* Variables should follow the same rules, but begin with a lower case letter: e.g. my\_variable\_name. etc. Variables should never be named something like aa, ab, ac, etc.
* We do not use CamelCase.
* Private members of classes should be postfixed with an underscore (e.g. var\_)
* Types should be postfixed with an \_t (e.g. My\_type\_t).
* Error messages should always list the name of the file generating the error and the line number from which the error message was printed. These can be set with the compiler Macros \_\_NAME\_\_ and \_\_LINE\_\_ so that they remain accurate as your source code changes.
* The GNU style guidelines recommend the use of gettext to make it easy to translate programs into different languages. We need to implement these standards in the future, however, we have not done so yet.
* Function definitions should have aligned braces:

Right:
	Void
	My_func (void)
	{
		return;
	}
Wrong:
	Void
	My_func (void) {
		return;
	}

* The return type should occur on the line immediately preceding the function name:

Right:
	int
	My_func (void)
	{
		return 1;
	}
Wrong:
	int My_func (void)
	{
		return 1;
	}

* Use the types defined in typedef.h when appropriate. For instance, when specifying a position in the genome, you should use the Id1\_t type.

Right:
	Id1\_t position;
Wrong:
	int position;

This prevents inconsistent usage for variables storing positions in the genome, which could be reasonably be declared as int, unsigned int, long int, etc. 

* 

##Coancestry Coefficients
##Maximum Likelihood
##Likelihoods
##Prior
##Posteriors
##Maximizing Posteriors
##Genotypic Correlation

Ultimately the genetic structure of a population is fully specified by the genotypes of the individuals that compose that population. This means that if we can accurately calculate all genotypic probabilities, then the calculation of any other population statistics becomes trivial. However, in order to calculate genotypic probabilities we must take account of the errors made in the genotyping process. Include inferring the presence of alleles that are not there, which can arise from sequencing error or mistakes made aligning to a reference, and failing to detect the presence of alleles which are there because of low coverage or biased sequencing of a single parental chromosome. We maximize a likelihood equation to account for sequencing error and the failure to sample genotypes, and then we test the data for fit to the parameters, and reject estimates where the data has a poor fit to the estimated parameters.

##Data types

##Tutorials
An introduction to quartets
Reading and writing to files
"Making likelihood functions"
"Maximizing a likelihood function"

