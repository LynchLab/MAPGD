
# MAPGD version 0.4.2

# Contents

[Introduction](https://github.com/LynchLab/src/mapgd_0.4/MAPGD#introduction)

[Basic Design](https://github.com/LynchLab/src/mapgd_0.4/MAPGD#basic-design)

[Style Guidelines](https://github.com/LynchLab/src/mapgd_0.4/MAPGD#style-guidelines)

[Tutorials](https://github.com/LynchLab/src/mapgd_0.4/MAPGD#tutorials)

[Task List](https://github.com/LynchLab/src/mapgd_0.4/MAPGD#task-list)

# Introduction

MAPGD is a series of related programs that estimate allele frequency, heterozygosity, Hardy-Weinberg disequilibrium and identity by descent (IBD) coefficients from population genomic data using statistically rigorous maximum-likelihood approach.

# Basic Design

* I/O
  * The user interface
   The user interacts with MAPGD through a [command line interface](https://en.wikipedia.org/wiki/Command-line_interface). One of the advantages of 
  * The SQL database
  * Data\_file
* Data
  MAPGD can most generally be described as a program which transforms data from one representation into a different representation. The features of the data which the user is interested in may not be obvious in the most basic representation of the sequencing data, and the user may wish to transform the data into a representation where the features are more obvious. Generally, ... 
  * Raw\_data
  * Data
* Indexing
  The standard format many programs use for referencing positions in the genome uses a string, used to reference a continuous region of the genome called either a [contig](https://en.wikipedia.org/wiki/Contig), scaffold or chromosome, and an integer, used to specify a position within that region. To reduce memory usage mapgd represent genome positions as a single integer number, called an absolute position. The File\_index class is used to convert between string/integer pairs and absolute positions. 

  * Indexed\_data
  * Indexed\_file
* Models
  In order to 
 
# Style Guidelines

mapgd is written to conform to the GNU style guidelines, at least to the extent that I have had time to read and implement the guidelines.

* Class names should be whole English words starting with an initial capital (to distinguish them from the core classes and types of c++): e.g My\_class\_name.
* Variables should follow the same rules, but begin with a lower case letter: e.g. my\_variable\_name. etc. Variables should never be named something like aa, ab, ac, etc.
* Do not to use CamelCase, where words  .
* Private or protected members of classes should be postfixed with an underscore (e.g. var\_)
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

# Models

In statistics, a [model](https://en.wikipedia.org/wiki/Statistical_model) is an idealized description of how the observed data was generated. This description includes uses ...

## Maximum Likelihood

## Likelihoods

## Priors

## Posteriors

## Maximizing Posteriors

# Genotypic Correlation

Ultimately the genetic structure of a population is fully specified by the genotypes of the individuals that compose that population. This means that if we can accurately calculate all genotypic probabilities, then the calculation of any other population statistics becomes trivial. However, in order to calculate genotypic probabilities we must take account of the errors made in the genotyping process. Include inferring the presence of alleles that are not there, which can arise from sequencing error or mistakes made aligning to a reference, and failing to detect the presence of alleles which are there because of low coverage or biased sequencing of a single parental chromosome. We maximize a likelihood equation to account for sequencing error and the failure to sample genotypes, and then we test the data for fit to the parameters, and reject estimates where the data has a poor fit to the estimated parameters.

# Tutorials

[An introduction to quartets]()
[Reading and writing to files]()
[Making models]()

# Task List

[ ] Write scripting tutorials for the main user page.

[ ] Claim mapgd.org.

[ ] Automate tutorial testing.

[ ] Automatic statistical and computational performance testing.

   Statistical test should report bais and MSE of MAPGD and several other programs (ANGSD/GATK/PLINK).
   Keep other programs up to date.
   Performance test should report ...

[ ] Write SQL tutorials.

[ ] Automatically generate figures from recent papers. 
