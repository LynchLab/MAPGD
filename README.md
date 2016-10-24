![Travis Says](https://travis-ci.org/LynchLab/MAPGD.svg?branch=master)
#[Download MAPGD](https://github.com/LynchLab/MAPGD/archive/master.zip)

##MAPGD version 0.4

[(C) Michael Lynch and Matthew Ackerman](https://github.com/LynchLab/MAPGD#-copyright-)

[<b>Visit the MAPGD development page</b>](https://lynchlab.github.io/MAPGD/)

###Contents 
####[Introduction](https://github.com/LynchLab/MAPGD#-introduction-)
#####[Quick start](https://github.com/LynchLab/MAPGD#-quick-start-)
#####[FAQ](https://github.com/LynchLab/MAPGD#-faq-)
#####[Changes](https://github.com/LynchLab/MAPGD#-changes-)
####Basic Usage
#####[Commands](https://github.com/LynchLab/MAPGD#-commands-)
#####[Common Options](https://github.com/LynchLab/MAPGD#-common-options-)
#####[Input/Output](https://github.com/LynchLab/MAPGD#-inputoutput-)
#####[Log-likelihood ratios](https://github.com/LynchLab/MAPGD#-log-likelihood-ratio-statistics-)
#####[Example Analysis](https://github.com/LynchLab/MAPGD#-example-analysis-)
#####[Other Useful Programs](https://github.com/LynchLab/MAPGD#-other-useful-programs-)
####Performance()
#####[Statistical](https://github.com/LynchLab/MAPGD#-statistical-performance-)
#####[Computational](https://github.com/LynchLab/MAPGD#-computational-performance-)
####Misc.
#####[IU users](https://github.com/LynchLab/MAPGD#-notes-for-indiana-university-users-)
#####[Sanger users](https://github.com/LynchLab/MAPGD#-notes-for-sanger-users-)
#####[References](https://github.com/LynchLab/MAPGD#-references-)

<h2> Introduction </h2>

MAPGD is a series of related programs that estimate allele frequency, heterozygosity, Hardy-Weinberg disequilibrium, linkage disequilibrium and identity-by-descent (IBD) coefficients from population genomic data using a statistically rigorous maximum likelihood approach. It is primarily useful for the analysis of low coverage population genomic data or for the analysis of pooled data (where many individuals are used to prepare a single sample).

<h2> Quick start </h2>

<h6>Installation</h6>

After clicking the "Download MAPGD.zip" button you will be prompted to save or open the file MAPGD-master.zip. Save this this file to the directory of your choice, then go to this directory in a terminal, for example /home/matthew/Downloads/

Then type:

	unzip MAPGD.zip
	cd MAPGD-master/
	make

The program can be installed for all users of a computer by typing:

	sudo make install

You will be prompted for your super-user password. If you do not have a super-user password for the system on which you are installing the software, you can type:

	make install DESTDIR=~/bin/

This will install the software in the ~/bin/ directory which should allow you to use the software by simply typing 'mapgd'. If not, add the following line to your .bashrc file:

	PATH=$PATH:~/bin

A quick test to make sure everything is working correctly can be conducted by typing:

	make test

This will output a of lines ending in PASS or FAIL to your terminal. Ideally all of the lines should say PASS. 

<h6>Mac installation</h6>

Mac users may not have developmental tools installed by default, or you may not have agreed to the xcode licence. You may have to install and configure xcode before using mapgd.
Once you have xcode you can type:

	make noomp.


<h6> Using mapgd </h6>

Mapgd works a number of [commands](https://github.com/LynchLab/MAPGD#-commands-) each with their own associated help.

 Generally you will start by creating mpileup files with samtools mpileup command, converting these mpileup files to the pro format with mapgd's proview command, and then beginning your analysis. Running the script "make\_and\_align\_reads.sh" in the src\test\ directory will simulate a genomics study and then run mapgd on the simulated data. You may also want to look at an [Example Analysis](https://github.com/LynchLab/MAPGD#-example-analysis-).

For analyzing individual labeled data you will probably want a command like this:

	samtools mpileup -q 25 -Q 25 -B individual1.sort.bam individual2.sort.bam 
	| mapgd proview -H seq1.header | mapgd allele 
	| mapgd filter -p 22 -E 0.01 -c 50 -C 200 -o allelefrequency-filtered.map

The filter command will limit the output to sites where the log-likelihood ratio of polymorphism is greater than 22 (-p 22), the error rate is less than 0.01 (-E 0.01) and the population coverage is between 50 and 200. 

And for analyzing pooled data you will probably want a command like this:

	samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam 
	| mapgd proview -H seq1.header | mapgd pool -a 22 -o allelefrequency-filtered.pol

The -a option will limit the output to sites where the log-likelihood ratio of polymorphism is greater than 22, which is a relatively stringent criteria.

In the case where the allele command is being used to estimated the seven genotypic correlation coefficients named pipes can be used for the I/O redirection.
	
	mkfifo map;
	mkfifo pro;
	samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam 
	| mapgd proview -H seq1.header | mapgd allele -p pro | mapgd filter -p 22 -E 0.01 -c 50 -C 200 > map;
	mapgd genotype -p pro -m map | mapgd relatedness > population-rel.out

And linkage disequilibrium can be calculated in a similar manner:

	mkfifo map;
	mkfifo pro;
	samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam 
	| mapgd proview -H seq1.header | mapgd allele -p pro | mapgd filter -p 22 -E 0.01 -c 50 -C 200 > map;
	mapgd linkage -p pro -m map > population-lnk.out

<h2> FAQ </h2>

<b> Why don't you provide information on indels? </b>

We currently do not have a likelihood model to account for errors in calling indels, so we cannot incorporate indels into our program at this time. This is on our TODO list, but may not occur for some time.

<b> How long does it take to run? </b>

Typical benchmarks with 16 threads on a 2.6 GHz put us at around 18,000 sites a second for 96 simulated individuals at 10x coverage. This means that the typical invertebrate population will take around three hours to analyze on a good computer, and a vertebrate genome might take a few days. If you have managed to sequence *Paris japonica* you're looking at three months of computation time if you run it on a single computer. However, mapgd is designed to be used in a cluster computing environment, and can make use of multiple nodes to dramatically reduce computation time. Running mapgd on 96 individuals with 150 Gbp genomes should be possible if a large number of nodes (say 50) are used. 

<b> Help, I can't get the program to compile. </b>

MAPGD requires a compiler that complies with the C++11 guidelines. If you have a C++11 compiler and mapgd will still not compile on your system, please e-mail me (Matthew Ackerman) so that I can work to correct the problem. I have compiled mapgd on:

* Ubuntu Linux, 14.04

* Red Hat Enterprise Linux 6.x

* A version of the Cray Linux Environment.

* OS X Yosemite.

To compile on OS X, you may need to type 'make noomp' because the default OS X does compiler does not support openmp. You may be able to obtain a compiler that supports openmp by typing:

	brew install gcc --without-multilib

* Windows 
mapgd is available on windows systems as a pair of binaries (mapgd-win32 and mapgd-win64) in the bin directory. These files are cross compiled with mingw, and are not extensively tested, so use at your own peril. The relatedness command is unavailable to windows users at the current time.  

<b> Help, the program keeps crashing/hanging </b>
 
The first thing you should do is e-mail me (for contact information type 'mapgd -h'). I will open a bug report and we can begin discussing how to fix the program. 

<b> How can I help? </b>

We have lots to do, and need plenty of help! 

* Create an Issue. 

	If you have any problems at all using this code, please click the !Issues button on the right hand menu to make us aware of the problem. Whether you can't compile and run the program, or locate one of the plethora of typos in the documentation and help files, please let us know! 

* Write a script to evaluate the statistical performance of a command. 

	One of the most important task we have right now is to compare the performance of mapgd to other programs that are available and make sure that our program is as good as it can be. Whether it is comparison of computational efficiency or statistical accuracy, we need a script to assess it! 

* Anything you can think of!

If you can do any of these things (or anything else that you think might help), you can contribute to the project by typing:

	git clone https://github.com/LynchLab/MAPGD/

This will create a clone of the repository that you can play around with on your own. Then, if you can make some change that might help you can mark the file to be changed by typing:

	git add FILENAME

Then type a short comment describing the change you have made using the command

	git commit -m "COMMENT"

Finally show me your changes by typing.

	git push

<h2> Changes </h2>

There have been a lot changes from 0.3. The format of input and output files has changed, and previous formats are no longer supported. The name of the 'ei' command has been changed to allele, and the 'ep' and 'cp' are now both part of the 'pooled' command. A standard file interface has been created (map-file) which handles all our reading and writing needs. The pro-file interface has been deprecated.  


<h2> Commands </h2>

Mapgd currently implements the following commands:

	allele                Estimates allele frequencies using individual data
	filter                Filter sites in '.map' files
	genotype              Calculate genotype probabilities for individuals
	linkage               Estimates linkage disequilibrium between loci
	pool                  Estimates allele frequencies using pooled data*
	proview               Prints data in the '.pro' file quartet format
	read                  Reads from an SQL database	
	relatedness           Estimates the 7 genotypic correlation coefficients 
	sam2idx               Reformats a sam header to an idx used by mapgd.
	vcf                   Converts output to the popular vcf format
	write                 Writes to an SQL database

Each command has a number of options that can be examined by the -h option. For example, to get a short help message you can type: 

	matthew@Ne:~$mapgd allele -h
	usage: mapgd allele  [--input] [--output] [--outpro] [--individuals] [--minerror] ... 

	mapgd allele version 0.4.1 written by Matthew Ackerman and Takahiro Maruki
	Uses a maximum likelihood approach to estimate population genomic statistics from an individually 'labeled' population.

	Options:
	  -i, --input		the input file for the program (default stdout).
	  -o, --output		the output file for the program (default stdin).
	  -p, --outpro		name of a 'cleaned' pro file (default none).
	  -I, --individuals	the individuals to be used in estimates.
				a comma seperated list containing no spaces, and the format X-Y can be used to specify a range (default ALL).
	  -m, --minerror	prior estimate of the error rate (defualt 0.001).
	  -H, --header		the name of a .idx file storing scaffold infomation
	  -M, --mincoverage	minimum coverage for an individual at a site for an individual to be used (default 4).
	  -g, --goodfit		cut-off value for the goodness of fit statistic (defaults 2.0).
	  -N, --number		cut-off value for number of bad individuals needed before a site is removed entirely (default 4).
	  -S, --skip		number of sites to skip before analysis begins (default 0).
	  -H, --noheader	disables printing a headerline.
	  -n, --newton		use newton-raphson likelihood maximization (slow but accurate).
	  -h, --help		prints this message
	  -v, --version		prints the program version

More detailed documentation for each command is being produced, and will be available shortly. 

<h2> Common Options </h2>

<b>-o	--output \<FILENAME\> </b> When this option is specified output file(s) with basename FILENAME will be created and with appropriate extensions will be added (possible file extensions are listed bellow). This option can be useful for generating files for human inspection because data will not be proceeded by a long list of scaffold names and lengths. However, several separate files will be need by some commands, which may make the management of files more difficult.

<b>-H	--header \<FILENAME\> </b> FILENAME will be treated as a header, and the .idx corresponding to basename will not be opened.

<b>-i   --input \<FILENAME\> </b> The file FILENAME will be opened for input, rather than taking input from stdin.

<h2> Input/Output </h2>

<b> Header lines </b> Ever file begins with two header lines, each beginning with the '@' character. The first header line list the name of the table in the SQL database in which data may be stored, along with the version of MAPGD used to create the table, and some formating information. The second header list the value stored in each column of the table. 

The table below lists some of the labels with their descriptions. It also lists the type of value stored in each column, but this will not be important for you unless you are writing a program which directly uses the binary output of mapgd. For a complete list of Labels and types see the file keys.txt in the source directory.

| Label	  | mapgd type 	| Description 							|
|:--------|:------------|:-------------------------------------------------------------:|
| MJ\_FREQ| float\_t	| frequency of the major allele					|
| MN\_FREQ| float\_t	| frequency of the minor allele					|
| MM\_FREQ| float\_t	| frequency of the major major genotype				|
| Mm\_FREQ| float\_t	| frequency of the major minor genotype				|
| mm\_FREQ| float\_t	| frequency of the minor minor genotype				|
| NULL\_ER| float\_t	| error rate assuming monomorphism				|
| ERROR   | float\_t	| maximum likelihood error rate					|
| HETERO  | float\_t	| heterozygosity of a site					|
| POLY\_LR| float\_t 	| log likelihood ratio of best fit/monomorphic			|
| HWE\_LR | float\_t 	| log likelihood ratio of best fit/Hardy–Weinberg  equilibrium	|
| SMPNUM  | size\_t	| sample number							|
| SMPNAME | std::string	| sample name							|
| SCFNAME | std::string	| the name of a scaffold, region of DNA				|
| POS     | id1\_t   	| position							|
| COVRAG  | count\_t 	| depth of coverage at a site					|
| IND\_INC| size\_t	| number of individuals used in a calculation			|
| IND\_CUT| size\_t	| number of individuals excluded from a calculation 		|
| EF\_CHRM| float\_t	| the effective number of chromosomes at a site			|
| GOF	  | float\_t	| goodness of fit value						|
| LENGTH  | id1\_t	| The length of a scaffold					|
| VERSION | std::string	| The version of mapgd used to make a file			|

Header lines can also contain an arbitrary (sanitized) string, which serves as a sample name, if appropriate. Below are example headers from each of the types of files produced by mapgd. 

<b> .idx files </b>

	@NAME:SCAFFOLDS	VERSION:0.4.1	FORMAT:TEXT
	@SCFNAME       	LENGTH

<b>.idx</b> files list the name and size of all the scaffolds in a reference genome. This file can be obtained from a .bam file using the samtools view -H command and reformatting the samtools header with the 'sam2idx' command. Idx files are automatically generated when running the proview command.

<b> .gof files </b>

	@NAME:SAMPLE	VERSION:0.4.1	FORMAT:TEXT
	@SMPNAME	GOF

<b>.gof</b> files are generated by the allele command. These short files list 'Goodness of fit' values for each sample in the population. These values can be used for filtering out samples that have been cross contaminated. 

<b> .map files </b> 

	@NAME:POSITIONS	VERSION:0.4.1	FORMAT:TEXT
	@SCFNAME    	POS	REF	MAJOR	MINOR	COVERAG	MJ_FREQ	MN_FREQ	ERROR	NULL_ER	F_STAT	MM_FREQ	Mm_FREQ	mm_FREQ	HETERO	POLY_LR	HWE_LR	GOF	EF_CHRM	IND_INC	IND_CUT	BEST_LL

<b>.map</b> files contain a list of estimated genotypic frequency obtained with the allele command. These files store test statistics for polymorphism and Hardy-Weinberg disequilibrium, as well as a small number of statistics which may prove useful for filtering variants, such as sequencing error rate and population depth of coverage.

<b> .pro files </b> 

	@NAME:QUARTETS	VERSION:0.4.1	FORMAT:TEXT
	@SCFNAME       	POS     REF     PA-001          PA-002          PA-003          ...
	scaffold_1      1       A       0/0/0/0         1/0/0/2         4/0/0/0

<b>.pro</b> files are the most basic input file for mapgd. These are plain text files containing three or more tab delimited columns. The first column is an arbitrary string which identifier a genomic region (e.g., a scaffold), the second column is an integer number specifying the location of a site on that scaffold, and the remaining column(s) contains four integer values separated by '/'s representing the number of times an A, C, G, and T was observed at the site (respectively).

We call this file format the .pro file format. Files in this format can be generated from mpileup files (that have been made without the -s and -O options by samtools mpileup) using the command "mapgd proview" 

If more than one bam file was used in the construction of the mpileup file, then these files will each appear as additional columns in the .pro File. The column names in mapgd default to the filename used to generate the column, and if more than one column is generated from a file, then the columns are numbered sequentially. Additionally, if multiple .pro or mpileup files are given as input to any command (proview included) these files can be merged for analysis (e.g "mapgd proview -i \*.mpileup" prints a merged pro file to the standard out). In order to preserve sample names for analysis, a ...

<b> .gcf files </b> 

The output of the genotype command. This stores the -log likelihood values that an individual is each of the three possible genotypes (Major Major, Major minor or minor minor) at each locus. 

<b> .rel files </b> 

The output of the relatedness command. This file stores the 7 genotypic correlation coefficients for all pairs of individuals and some log likelihood ratio test statistics. 

<b> .pol files </b> 

	@NAME:SAMPLE	VERSION:0.4.1	FORMAT:TEXT
	@SCFNAME       	POS     MAJOR   MINOR   COVRAG  ERROR   Sample_0        Sample_1        Sample_2        Sample_3        Sample_4        Sample_5        Sample_6        Sample_7
	scaffold_3      1       T       A       3       0.001   .../.../.../... 0.6/0.3/24./28.	.../.../.../...	.../.../.../...	.../.../.../...	.../.../.../...	.../.../.../...	.../.../.../...	

The output of the pooled command. This stores best estimates of allele frequencies and log likelihood ratio test for polymorphism and fixed substitutions at each locus. The log likelihood ratios are Polymorphic/Fixed Major, Polymorphic/Fixed Minor, Fixed Major/Fixed Minor. Because polymorphism has one more free parameter than the fixed states, log likelihood ratios for polymorphism will always be positive. For the Fixed Major vs. Fixed Minor, the statistic can be either positive (in which case it is more likely that the sample is fixed for the Major allele) or negative (in which case it is more likely that the sample is fixed for the minor allele. 

So, for example, in Sample\_1, the maximum likelihood estimate of the allele frequency is 0.6, but the log likelihood ratio test against against the fixed major allele is only 0.3, which is not significant. 
<h3> The SQL database </h3>

MAPGD includes the ability to direct all output to an SQL database. This decrease storage space by eliminating redundant information from files (such as the initial position label in all indexed files) as well as decreasing the number of separate files saved to disk. Currently the tables are:

| TABLES 	| Data 	 | Primary Keys |
|:--------------|:-------|:------------:|
| SCAFFOLDS	|	 |
| POSITIONS	|	 |			|
| SAMPLE	|	 |			|
| REGIONS	|	 |Scaffold, Start, Stop	|
| SAMPLE\_PAIRS	|

In order to 

<h2> Log-likelihood ratio statistics </h2>

Most of the commands in mapgd report log-likelihood ratio statistics. These statistics should be chi-square distributed. The number of degrees of freedom of the statistic depend on the number of parameters being estimated. In the case of pooled population data there is one degree of freedom, for the allele polymorphic statistic there are two, and for the Hardy-Weinberg equilibrium statistic there is one. For the relatedness statistics there is one degree of freedom for each parameter, and seven degrees of freedom between the best-fit and null statistic. Significance at the 0.05, 0.01, 0.001 levels requires that the likelihood-ratio test statistic exceed 3.841, 6.635, and 10.827, respectively. Please consider including a correction for multiple testing if you wish to limit the number of type I errors in your data set. Other critical values can be obtained in R by typing : chisq(VALUE, df=DEGREES OF FREEDOM).  

<h2> Extras </h2>

<h6> pedigree_calc.py </h6>A simple script to predict the coefficients of IBD from merlin formated pedigree files.
<h6> simulate_genomic_study.sh </h6> A script that generates raw reads from a population and a reference file. This script was used to compare the performance of mapgd with similar programs. 
<h6> make_bordyen_theta.py  </h6> A script to automatically generate some source code for minimization. The code does not work correctly. 
<h6> make_newton_rho.py </h6>  A script to automatically generate some source code for minimization. The code does not work correctly.
<h6> make_newton_theta.py </h6> A script to automatically generate some source code for minimization. The code doe not work correctly.
<h6> marker.txt </h6> A list of 30,000 some odd markers for asexuality that were found in Tucker et al.
<h6> score_markers.py </h6> A script to record the % of makers in a file that and individual has. Usage python score\_markers.py marker.txt GENOTYPES.gcf.

<h2> Example Analysis </h2>

To begin any of these analyses, a ".pro" file must be created from a mpileup file. This is done with the proview command, which needs one or more mpileup files and a single index file.

For example, if the sequencing center gives you two files called "seq1.fastq" and "seq2.fastq" your entire work flow might look something like this: 

	bwa aln Reference.fna seq1.fastq > seq1.sai
	bwa sam Reference.fna seq1.sai seq1.fastq > seq1.sam
	samtools view -bS seq1.sam > seq1.bam
	samtools sort seq1.bam seq1.sort
	samtools index seq1.sort.bam

You would then type some similar commands to map reads from seq2 to the same reference.

Next the data have to be converted into some format the mapgd can read. This is done by creating an 'mpileup' file with samtools:

	samtools mpileup -q 25 -Q 25 -B seq1.sort.bam seq2.sort.bam > population.mpileup

getting the samfile header by typing:

	samtools view -H seq1.sort.bam > seq1.header

and then converted the mpileup file to a .pro file: 

	mapgd proview -i metapopulation.pileup -H seq1.header > metapopulation.pro

If the .pro file contains pooled (e.g. individuals cannot be distinguished) data you will want to run the pooled command:  

	mapgd pool -i metapopulation.pro -o population1_allelfrequencies.map

Alternatively, if the .pro file contains individual data you will want to run the allele command:

	mapgd allele -i metapopulation.pro -o population1_allelfrequencies.map

One advantage of the above work flow is that each of the files can be inspect visually to explore the data. If you do not need to do this then it will be faster to create binary data and use I/O redirection.

<h5> The above work flow using  I/O redirection </h5>

For analyzing individual labeled data you will probably want a command like this:

	samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam 
	| mapgd proview -H seq1.header | mapgd allele | mapgd filter -p 22 -E 0.01 -c 50 -C 200 -o allelefrequency-filtered.map

And for analyzing pooled data you will probably want a command like this:

	samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam 
	| mapgd proview -H seq1.header | mapgd pool -a 22 -o allelefrequency-filtered.pol

In the case where the allele command is being used to estimated the seven genotypic correlation coefficients named pipes can be used for the I/O redirection.
	
	mkfifo map;
	mkfifo pro;
	samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam 
	| mapgd proview -H seq1.header | tee pro | mapgd allele | mapgd filter -p 22 -E 0.01 -c 50 -C 200 > map;
	mapgd genotype -p pro -m map | relatedness.py -o population.rel
<h3> Statistical Performance </h3>

![Figure1](https://github.com/LynchLab/MAPGD/extras/automated_figures/Ackerman2016b/figure1.jpg)

<h3> Computational Performance </h3>

![Figure2](https://github.com/LynchLab/MAPGD/extras/automated_figures/Ackerman2016b/figure2.jpg)

<h3> Other Useful Programs </h3>

This program is intended for use with ".pro" described in the previous section. A slightly modified version of the program sam2pro, written by Bernhard Haubold, is included in this package and can be run by typing "mapgd proview". Programs for converting other file formats to ".pro" files are available in stand alone form http://guanine.evolbio.mpg.de/mlRho/

The workflow described above requires the programs bwa and samtools.

To download bwa please visit http://bio-bwa.sourceforge.net/

To download samtools please visit http://www.htslib.org/

<h5> For windows users </h5>

Windows is currently unsupported, but you may try to compile the code and fix it yourself. I have tried to refrain from using any platform specific libraries, so it may not be too much work.

<h2> Notes for any High Performance Computer users</h2>
MAPGD is able to take advantage of mutli-threading and cluster computing environment. Make sure you know how to submit jobs that execute on multiple nodes with mutliple CPUs. To test whether your jobs are running correctly you can run the script 'speedtest' in the src/test directory. Ideally you should roughly linear gains from adding threads and CPUs up to the maximum avalible to you. 

<h2> Notes for Indiana University users </h2>
When submitting PBS scripts please *make sure to specify the number of threads to use with the ppn option* 

	#PBS -l ppn=16

*Bigred2* has several different programming environments. To compile the code you will have to type 

	module rm PrgEnv-cray
	module load PrgEnv-gnu 
	module load gsl/1.15	

*Karst* may require loading

	module load gsl/1.15	

*Mason* will require 

	module rm gcc
	module load gcc/4.9.2
	module load gsl/1.15	

<h2> Notes for Sanger users </h2>

*Farm3* will require

	configure 'CXX=/software/gcc-4.9.2/bin' 'CXXFLAGS=-static'
	make	

Word of warnding, the static linking is likely to make the program run slower, but I'm not sure by how much. If you are able you may want to change your environmental variables so that 'configure; make' works. 

<h3> References </h3>

Ackerman, M. S., T. Maruki and M. Lynch. "MAPGD a program for the maximum likelihood analysis of population data." In prep.

For output of the allele command please cite:

Maruki, T., and M. Lynch. 2015  "Genotype-Frequency Estimation from High-Throughput Sequencing Data." Genetics 201.2: 473-486.

For output of the pool command please cite:

Lynch, M., D. Bost, S. Wilson, T. Maruki, and S. Harrison. 2014  "Population-Genetic Inference from Pooled-Sequencing Data." Genome Biol Evol 6:1210-1218.

For the output of the linkage command please cite:
	
Maruki, T., and M. Lynch 2014 "Genome-Wide Estimation of Linkage Disequilibrium from Population-Level High-Throughput Sequencing Data." Genetics 197: 1303-1313;

For output of the relatedness command please cite:

Ackerman, M. S., P. Johri, K. Spitze and M. Lynch, 2015  A general statistical model for coefficients of relatedness and its application to the analysis of population-genomic data. In prep. 

<h3> Copyright </h3>

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

For a copy of the GNU General Public License write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
