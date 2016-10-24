
# MAPGD version 0.4.2

# Contents

[Introduction](https://github.com/LynchLab/src/MAPGD_0.4/MAPGD#introduction)

[Basic Design](https://github.com/LynchLab/src/MAPGD_0.4/MAPGD#basic-design)

[Style Guidelines](https://github.com/LynchLab/src/MAPGD_0.4/MAPGD#style-guidelines)

[Tutorials](https://github.com/LynchLab/src/MAPGD_0.4/MAPGD#tutorials)

[Task List](https://github.com/LynchLab/src/MAPGD_0.4/MAPGD#task-list)

# Introduction

MAPGD is a series of related programs that estimate allele frequency, heterozygosity, Hardy-Weinberg disequilibrium and identity by descent (IBD) coefficients from population genomic data using statistically rigorous maximum-likelihood approach.

# Basic Design
* The user interface
  The user interacts with MAPGD through a [command line interface](https://en.wikipedia.org/wiki/Command-line_interface). One of the advantages of a CLI is that it allows users to easily document and automate their analysis, however a CLI often provides less contextual help to users and users may be generally unfamiliar with CLI. We hope to address this by enforcing rich feature documentation through the use of MAPGD's [interface](https://lynchlab.github.io/MAPGD/interface_8h_source.html), which requires us to prove short verbal descriptions of all options and formats contextual help. MAPGD also provides documentation through man pages. 

* Data
  MAPGD can most generally be described as a program which transforms data from one representation into a different representation. The features of the data which the user is interested in may not be obvious in the most basic representation of the sequencing data, and the user may wish to transform the data into a representation where these features are more obvious. We attempt to organize data into minimal classes that associate one particular representation . . .  This choice 
  * Raw\_data
  * Data
* Indexing
  The standard format many programs use for referencing positions in the genome uses a string, used to reference a continuous region of the genome called either a [contig](https://en.wikipedia.org/wiki/Contig), scaffold or chromosome, and an integer, used to specify a position within that region. To reduce memory usage MAPGD represent genome positions as a single integer number, called an absolute position. The File\_index class is used to convert between string/integer pairs and absolute positions. 
* File interface
  MAPGD provides a basic file interface for reading and writing data. This interface allows the user to 
  * Indexed\_data 
  Data which is logically associated with a position in the genome can be declared as Indexed\_data. This allows 
  * Doubly\_indexed\_data 
* I/O
  MAPGD enforces a few basic rules with respect to IO. The input of most commands should be a File, and the output should be a File. Commands that naturally take more than one Data\_type as input should define a single Data\_type that unites these Data\_types so that they can be written to a single file. Then, a command that specifically takes the different Data\_types as input and outputs the single new Data\_type should also be implemented. This ensures that MAPGD can effectively use the SQL database back end.  
  * The SQL database
* Models
  In order to 
 
# Style Guidelines

MAPGD is written to conform to the GNU style guidelines, at least to the extent that I have had time to read and implement the guidelines.

* Class names should be whole English words starting with an initial capital (to distinguish them from the core classes and types of c++): e.g My\_class\_name.
* Variables should generally follow the same rules, but begin with a lower case letter: e.g. my\_variable\_name. etc. Variables should never be named something like aa, ab, ac, etc.
* Common abbreviations are acceptable for variable names (but not for class or function names). For instance, tmp would be preferred over temporary, etc. Use your own judgment on whether the use of a variable is obvious enough from context that a long descriptive name is unnecessary (loop counters, temporary storage, etc.). When in  doubt, error on the side of excess description. 
* Do not to use [CamelCase](https://en.wikipedia.org/wiki/Camel_case), where several words are joined together and each word begins with a capital letter. Instead join words with the underscore, since non-native English speakers allegedly find this easier to read.
* Private or protected members of classes should be postfixed with an underscore (e.g. var\_). 
* Types should be postfixed with an \_t (e.g. My\_type\_t).
* Avoid type prefixes.
  Prefixes which encode type information (e.g. ui\_count, to show that count is an unsigned int) should generally not be used. In most IDEs it is quite easy to display  types, and at any rate, if the code is reasonably well organized it shouldn't be too difficult to find the declaration. 
* Error messages should always list the name of the file generating the error and the line number from which the error message was printed. These can be set with the compiler Macros \_\_NAME\_\_ and \_\_LINE\_\_ so that they remain accurate as your source code changes.
* The GNU style guidelines recommend the use of gettext to make it easy to translate programs into different languages. We need to implement these standards in the future, however, we have not done so yet. 
* Function definitions should have aligned braces:

   Right:
    void
    My_func (void)
    {
       return;
    }
   Wrong:
    void
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


* Do not assume that characters will be encoded in ASCII. We would like MAPGD to work for international users, who may have text encoded in UTF-8 etc. Please do not break this.
* Keep source file small. While you do not have to ensure that only a single object exists in source file style, they should certainly not contain logically independent objects. This allows the code to easily be reorganized by moving files around, allows for re-use of objects in other programs, and generally lowers the barrier for other programmers to contribute  to the project.
* Each object files should compile independently. 
* Comment your code. 
  We will be using [DOXYGEN](http://www.stack.nl/~dimitri/doxygen/) to comment the code. Header files should include descriptions for essentially every line of code that will appear in the development website. Comments may be more sparse in the source file, since other developers have less need to poke around in these files as long as everything is working properly. Still, write your code to be read by other people.  

# Models

In statistics, a [model](https://en.wikipedia.org/wiki/Statistical_model) is an idealized description of how data is generated. If a statistical model of data exist, it can be used to estimate parameters used to generate the data through the process of likelihood maximization.

## Likelihoods

## Maximum Likelihood

## Priors

## Posteriors

## Maximizing Posteriors

## Detecting violations of the model

# Useful representations of data

## Allele frequencies

## Genotypic likelihoods

## Genotypic Correlation

Ultimately the genetic structure of a population is fully specified by the genotypes of the individuals that compose that population. This means that if we can accurately calculate all genotypic probabilities, then the calculation of any other population statistics becomes trivial. However, in order to calculate genotypic probabilities we must take account of the errors made in the genotyping process. Include inferring the presence of alleles that are not there, which can arise from sequencing error or mistakes made aligning to a reference, and failing to detect the presence of alleles which are there because of low coverage or biased sequencing of a single parental chromosome. We maximize a likelihood equation to account for sequencing error and the failure to sample genotypes, and then we test the data for fit to the parameters, and reject estimates where the data has a poor fit to the estimated parameters.

### Jaquard's condensed coefficients of identity

## Zygosity Correlations

The genotypic correlation between two loci within a population is generally described in terms of linkage disequilibrium, but can more generally be seen as second order zygosity correlations between loci. These correlations exist within individuals, and ultimately linkage disequilibrium is simply the average second order zygosity correlation across individuals within the population.

# Tutorials

[An introduction to quartets]()
[Reading and writing to files]()
[Making models]()

# Task List

[x] Claim MAPGD.org.

   The website should probably just redirect to the github repository.

[ ] Implement gettext for internationalization.

[ ] Get SQL back end working again.

[ ] Get MPI working again.

[ ] Write scripting tutorials for the main user page.

   Write a tutorial to show 

[ ] Write SQL tutorials.

   SQL tutorials should show users how to quickly calculate various summary statistics form the SQL database, such as piN and piS.

[ ] Automate tutorial testing.
   
   Both coding tutorials and example work flows in the README should be automatically tested before each push to the main repository.

[ ] Automatic statistical and computational performance testing.

   Statistical test should report bais and MSE of MAPGD and several other programs (ANGSD/GATK/PLINK).
   Add scripts to keep other programs up to date.
   Performance test should look at how computations times scale with input size, number of threads, and number of nodes.
   Results of Statistical and Performance test should automatically be displayed in the readme.


[x] Automatically generate figures from recent papers. 

   These figures and scripts are placed in the directories Ackerman 2016a and Ackerman 2016b.
