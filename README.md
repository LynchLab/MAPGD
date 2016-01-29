##This readme is being update and does not currently accurately describe the program.

##MAPGD version 0.4

[<b>Download MAPGD.zip</b>](https://github.com/LynchLab/MAPGD/archive/master.zip) [(C) Michael Lynch](https://github.com/LynchLab/MAPGD#-copyright-)

[<b>Visit the MAPGD development page</b>](https://lynchlab.github.io/MAPGD/)

###Contents 
####[Introduction](https://github.com/LynchLab/MAPGD#-introduction-)
#####[FAQ](https://github.com/LynchLab/MAPGD#-faq-)
#####[Changes](https://github.com/LynchLab/MAPGD#-changes-)
#####[Installation](https://github.com/LynchLab/MAPGD#-installation-)
####Basic Usage
#####[Commands](https://github.com/LynchLab/MAPGD#-commands-)
#####[Input/Output](https://github.com/LynchLab/MAPGD#-inputoutput-)
#####[Examples](https://github.com/LynchLab/MAPGD#-examples-)
####Misc.
#####[Extras](https://github.com/LynchLab/MAPGD#-extras-)
#####[Bigred2](https://github.com/LynchLab/MAPGD#-bigred2-)
#####[Mason](https://github.com/LynchLab/MAPGD#-Mason-)
#####[Karst](https://github.com/LynchLab/MAPGD#-karst-)
#####[References](https://github.com/LynchLab/MAPGD#-references-)

<h3> Introduction </h3>

MAPGD is a series of related programs that estimate allele frequency, heterozygosity, Hardy-Weinberg disequilibrium, linkage disequilibrium and identity by descent (IBD) coefficients from population genomic data using statistically rigorous maximum likelihood approach. It is primarily useful for the analysis of low coverage population genomic data or for the analysis of pooled data (where many individuals are used to prepare a single sample). Although other tools may give similar statistical estimates with less computational investment when coverage is high, the IBD estimates of other programs may not converge even at high coverage.

<h5> FAQ </h5>

<b> Why don't you provide information on indels? </b>

We currently do not have a likelihood model to account for errors in calling indels, so we cannot incorporate indels into our program at this time. This is on our TODO list, but may not occur for some time.

<b> How long does it take to run? </b>

Typical benchmarks with 16 threads on a 2.6 GHz put us at around 18,000 sites a second for 96 simulated individuals at 10x coverage. This means that the typical invertebrate population will take around three hours to analyze on a good computer, and a vertebrate genome might take a few days. If you have managed to sequence *Paris japonica* you're looking at three months of computation time. We are not currently focused on scaling our program up for larger genomes, however, if you want to analyze something like *japonica*, let us know and we will make MPS and CUDA priorities! 

<b> Help, I can't get the program to compile. </b>

MAPGD has requires a compiler that complies with the C++11 guidelines, however, I am happy to work with you to ensure that MAPGD can run on a wide variety of systems. I have compiled MAPGD on:

* Ubuntu Linux, 14.04

* Red Hat Enterprise Linux 6.x

* A version of the Cray Linux Environment.

* OS X Yosemite.

To compile on OS X, you need to type 'make noomp' because the default OS X does compiler does not support openmp.  
 
Version 0.1 will compile on windows 32 and 64 bit versions, but getting current versions of MAPGD working on windows is near the bottom of the TODO list.

Long story short: If you have problems, please please please e-mail me, and I will work to get it to compile for you as quickly as I can. 

<b> Help, the program keeps crashing/hanging </b>
 
The first thing you should do is e-mail me (for contact infromation type 'mapgd -h'). I will open a bug report and we can begin discussing how to fix the program. 

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

<h5> Changes </h5>

There have been a lot changes from 0.3. The format of input and output files has changed, and previous formats are no longer supported. The name of the 'ei' command has been changed to allele, and the 'ep' and 'cp' are now both part of the 'pooled' command. A standard file interface has been created (map-file) which handels all our reading and writing needs. The pro-file interface has been depricated.  
 
<h5> Commands </h5>

Mapgd currently implements the following commands:
* allele                Estimates allele frequencies using individual data
* filter                Filter sites in '.map' files
* genotype              Calculate genotype probabilities for individuals
* linkage               Estimates linkage disequilibrium between loci
* pooled                Estimates allele frequencies using pooled data
* proview               Prints data in the '.pro' file quartet format
* relatedness           Estimates the 7 IBD coefficients 
* vcf                   Converts output to the popular vcf format

In the near future we hope to implement the commands:
* read			Reads from an SQL database	
* write			Writes to an SQL database

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
				a comma seperated list containing no spaces, and the format X-Y can be used to specify a range (defualt ALL).
	  -m, --minerror	prior estimate of the error rate (defualt 0.001).
	  -H, --header		the name of a .idx file storing scaffold infomation
	  -M, --mincoverage	minimum coverage for an individual at a site for an individual to be used (defualt 4).
	  -g, --goodfit		cut-off value for the goodness of fit statistic (defaults 2.0).
	  -N, --number		cut-off value for number of bad individuals needed before a site is removed entirely (default 4).
	  -S, --skip		number of sites to skip before analysis begins (default 0).
	  -H, --noheader	disables printing a headerline.
	  -n, --newton		use newton-raphson likelihood maximization (slow but accurate).
	  -h, --help		prints this message
	  -v, --version		prints the program version

More detailed documentation for each command is being produced, and will be availible shortly. 

<h5> Input/Output </h5>

<b> Header lines </b> Ever file begins with two header lines, each begining with the '@' character. The first header line list the name of the table in the SQL database, along with the version of MAPGD used to create the table. The second header list the value stored in each column of the table. 

The Table bellow lists potential column labels with their discriptision. It also lists the type of value stored in each column. This information is provided for developers, but will not be important unless you are writing a program which directly uses the binary output of mapgd.

| Label	  | mapgd type 	| Description 							|
|:--------|:------------|:-------------------------------------------------------------:|
| KEY  	  | char[7] 	| a unique label for data stored in columns			|
| TYPE    | char[7] 	| the type of data stored in KEY columns			|
| DESC 	  | std::string	| string  a verbal description of data stored by the key	|
| GENOPRB | Genotype	| A set of three genotypic probabilities			|
| GENOTYP | gt\_t  	| genotype code							|
| MJ\_FREQ| float\_t	| frequency of the major allele					|
| MN\_FREQ| float\_t	| frequency of the minor allele					|
| MM\_FREQ| float\_t	| frequency of the major major genotype				|
| Mm\_FREQ| float\_t	| frequency of the major minor genotype				|
| mm\_FREQ| float\_t	| frequency of the minor minor genotype				|
| NULL\_ER| float\_t	| error rate assuming monomorphism				|
| ERRROR  | float\_t	| maximum likelihood error rate					|
| HETERO  | float\_t	| heterozygosity of a site					|
| POLY\_LR| float\_t 	| log likelihood ratio of best fit/monomorphic			|
| HWE\_LR | float\_t 	| log likelihood ratio of best fit/Hardy Weinburg eqaulibrium	|
| SMPNUM  | size\_t	| sample number							|
| SMPNAME | size\_t	| sample name							|
| SCFNAME | id0\_t 	| the name of a scaffold, a ? numbered region of DNA.		|
| POS     | id1\_t   	| position							|
| COVRAG  | count\_t 	| depth of coverage at a site					|
| IND\_TOT| size\_t 	| number of samples with data					|
| IND\_INC| size\_t	| number of individuals used in a calculation			|
| IND\_CUT| size\_t	| number of individuals excluded from a caclculation 		|
| EF\_CHRM| float\_t	| the effective number of chromosomes at a site			|
| SOMATIC | bool	| a flag to indicate that the record is a somatic mutation	|
| VALID   | bool	| a flag to indicate that a record has been validated?		|
| STRNDBS | float\_t	| strand bias at this position					|
| MAPQ0	  | count\_t	| number of MAPQ==0 reads covering this position		|
| HAPMAP2 | bool	| membership in hapmap2						|
| END     | size\_t	| ending row of a variant described by this record		|
| ANCSTRL | gt\_t	| ancestral allele						|
| GOF	  | float\_t	| goodness of fit value						|
| DEFAULT | bool	| a flag to indicate that a record represents a default value	|
| QUART   | Quartet	| a set of counts of nucleotides at a position 			|
| LENGTH  | id1\_t	| The length of a scaffold					|
| VERSION | std::string	| The version of mapgd used to make a file			|

Header lines can also contain an arbitrary (sanitized) string, which serves as a sample name, if appropriate. Below are example headers from each of the types of files produced by mapgd. 

<b> .idx files </b>

	@NAME:SCAFFOLDS	VERSION:0.4.1
	@SCFNAME       	LENGTH

file that list the name and size of all the scaffolds in a referance genome. This file can be obtained from a .bam file using the samtools view -H command, and reformating the samtools header with the 'sam2idx' command. Idx files are automatically generated when running the proview command.
Index files store the length (LENGTH) of each scaffold (SCFNAME) in the genome. This information is stored in a table named SCAFFOLDS in the database.

<b> .gof files </b>

	@NAME:SAMPLE	VERSION:0.4.1
	@SMPNAME	GOF

Generated by the allele command. This short file list 'Goodness of fit' values for each sample in the population. These values can be used for filtering out samlpes that have been cross contaminated. 
.gof files store the goodness of fit score (GOF) for each sample (SMPNAME) in the population. This information is stored in a table named SAMPLE in the database.

<b> .map files </b> 

	@NAME:POSITIONS	VERSION:0.4.1
	@SCFNAME    	POS	REF	MAJOR	MINOR	COVERAG	MJ\_FREQ	MN\_FREQ	ERROR	NULL\_ER	F\_STAT	MM\_FREQ	Mm\_FREQ	mm\_FREQ	HETERO	POLY\_LR	HWE\_LR	GOF	EF\_CHRM	IND\_INC	IND\_CUT	BEST\_LL

The output for genotypic frequency estiamtes obtained with the allele command. This file storages test statistics for polymorphism and Hardy-Weinberg disequilibrium, as well as a small number of statistics which may prove useful for filtering variants, such as sequencing error rate and population depth of coverage.

.pro files store the goodness of fit score (GOF) for each sample (SMPNAME) in the population. This information is stored in a table named SAMPLE in the database.

| TABLES 	|
|:--------------|:-----:|
| SCAFFOLDS	|
| POSITIONS	|
| SAMPLE	|
| REGIONS	|
| SAMPLE\_PAIRS	|

<b> .pro files </b> 

	@NAME:QUARTETS	VERSION:0.4.1
	@SCFNAME       	POS	REF	PA-001		PA-002		PA-003		...
	scaffold\_1	1	A	000/000/000/000	001/000/000/002	004/000/000/000

The most basic input file is a plain text file consists of three tab delimited columns : The first column is an arbitrary string which identifier a genomic region (e.g., a scaffold), the second column is an integer number specifying the location of a site on that scaffold, and the final column contains four integer values separated by '/'s representing the number of times an A, C, G, and T was observed at the site (respectively). The string-integer pair must be unique.  

We call this file format the .pro file format. Files in this format can be generated from mpileup files (that have been made without the -s and -O options by samtools mpileup) using the command "mapgd proview" 

If more than one bam file was used in the construction of the mpileup file, then these files will each appear as additional columns in the .pro File. The column names in mapgd default to the filename used to generate the column, and if more than one column is generated from a file, then the columns are numbered sequenctially. Additionally, if multiple .pro or mpileup files are given as input to any command (proview included) these files can be merged for analysis (e.g “mapgd proview -i \*.mpileup” prints a merged pro file to the standard out).

<b> .gcf files </b> The output of the genotype command. This stores the -log likelihood values that an individual is each of the three possible genotypes (Major Major, Major minor or minor minor) at each locus. 

<b> .rel files </b> The output of the relatedness command. This file stores the 7 genotypic correlation coefficents for all pairs of individuals and some log likelihood ratio test statistics. 


Under the assumption of a chi-square distribution for the test statistic with one degree of freedom, significance at the 0.05, 0.01, 0.001 levels requires that the likelihood-ratio test statistic exceed 3.841, 6.635, and 10.827, respectively. 

In principle, the 95% support interval can be obtained by determining the changes in the estimate of the minor allele frequency in both directions required to reduce the log likelihood by the appropriate chi-square value (e.g., 3.841) although this is not currently implemented. 

By default the program prints information to the file "dataout.txt" and this file will appear in the same location as the program. If an alternative file name is desired simply type if "mapgd -o FILENAME" where FILENAME is the name of your output file.

<b> The output of cp: </b> Columns 1 and 2 are site identifiers (ID1 and ID2); 3 and 4 designate major and minor nucleotides (major and minor);
For each population two columns are printed : The maximum likelihood estimate of the major allele frequency in that population (freq\_P), and the test statistic whether this frequency differes from the average frequency across all samples;

The final two columns are the major allele frequency in the metapopulation (meta\_P) and the maximum likelihood estimate of the error rate (Error). Output columns are tab delimited.

By default the program prints information to the file "dataout.txt" and this file will appear in the same location as the program. If an alternative file name is desired simply type if "mapgd -o FILENAME" where FILENAME is the name of your output file.

<b> The output of ei: </b> Columns 1 and 2 are site identifiers (ID1 and ID2); 3 and 4 designate major and minor nucleotides (major and minor); 
For each population in the sample four columns are printed : a major allele frequency (freq\_P), test statistic for polymorphism (ll\_poly), test statistic for fixed for the minor allele (ll\_fixed), and coverage;
The final coloumn is the maximum likelihood estimate of the error rate (Error).

Under the assumption of a chi-square distribution for the test statistic with one degree of freedom, significance at the 0.05, 0.01, 0.001 levels requires that the likelihood-ratio test statistic exceed 3.841, 6.635, and 10.827, respectively. 

In principle, the 95% support interval can be obtained by determining the changes in the estimate of the minor allele frequency in both directions required to reduce the log likelihood by the appropriate chi-square value (e.g., 3.841) although this is not currently implemented. 

By default the program prints information to the file "dataout.txt" and this file will appear in the same location as the program. If an alternative file name is desired simply type if "mapgd -o FILENAME" where FILENAME is the name of your output file.

<b> .gcf files </b> Again, this format is actively being developed, but .gcf files are intended to serve as a close analog of vcf files, with the only intentional difference being the storage of genotpyic likelihoods as floating point numbers.

<h4> Extras </h4>

<h6> mapgdutils.py </h6>The script mapgdutils.py is a stop gap measure used to generate genotypic likelihoods (i.e. gcf files) for the analysis of the relatedness of a pair of individuals. A script of this name is likely to persist in future version of the program to perform computationally simple task that may be useful; however, we hope to generation genotypic likelihoods directly in mapgd in the near future.

	mapgd proview -i \*.mpileups | mapgd allele -p FILENAME.pro -o FILENAME.map 
	python mapgdutils.py --pro FILENAME.pro --map FILENAME.map > GENOTYPES.gcf

<h6> pedigreecalc.py </h6>A simple script to predict the coefficients of IBD from merlin formated pedigree files.

<h6> relatedness.py </h6> A script to esimate the IBD coefficients of pairs of individuals from the '.gcf' files generated by mapgduils.py. To run relatedness.py type python relatedness.py 

<h6> simulatempileup.py </h6> A script that generates raw reads from a population and a reference file. This script was used to compare the performance of mapgd with similar programs. 

<h6> .py </h6> A script that generates raw reads from a population and a reference file. This script was used to compare the performance of mapgd with similar programs. 
<h6> MakeBordyenTheta.py  </h6> A script to automatically generate some source code for minimization. 
<h6> MakeNewtonRho.py </h6>  A script to automatically generate some source code for minimization. 
<h6> MakeNewtonTheta.py </h6> A script to automatically generate some source code for minimization. 
<h6> marker.txt </h6> A list of 30,000 some odd markers for asexuallity that were found in Tucker et al.
<h6> score_markers.py </h6> A script to record the % of makers in a file that and individual has. Usage python score\_markers.py marker.txt GENOTYPES.gcf.

<h5> Making the Input File </h5>

Currently mapgd allows for the estimation of allele frequencies in individually labeled and pooled population sequence through the "allele" and "pooled" command respectively, the comparison of allele frequencies between two populations with the "cp" command,  the conversion of mpileups to the ".pro" format though the proview command, and the estimatation of seven components of pair-wise relatedness through the "rel" command.

For example, if the sequencing center gives you two files called "seq1.fastq" and "seq2.fastq" your entire work flow might look something like this: 

<h4> Examples </h4>

<h6> Getting Fst for two pooled populations </h6> Generating Fst estimates for popualtions sequenced with pooled sequenc

First you will have to map your reads to some reference. This could be done with the program bwa by typing:
 
	bwa aln Reference.fna seq1.fastq > seq1.sai
	bwa sam Reference.fna seq1.sai seq1.fastq > seq1.sam
	samtools view -bS seq1.sam > seq1.bam
	samtools sort seq1.bam seq1.sort
	samtools index seq1.sort.bam

You would then type some similar commands to map reads from seq2 to the same reference.

Next the data have to be converted into some format the mapgd can read. Currently the easiest route to accomplish this is to create an 'mpileup' file with samtools:

	samtools mpileup -q 25 -Q 25 -B seq.sort.bam seq.sort.bam > population.mpileup

Then mpileup file must be converted to a .pro file which mapgd can use directly: 

	mapgd proview -i metapopulation.pileup > metapopulation.pro


Finally analysis can be run on the data. If the .pro file contains population data you might run either the ep command or the cp command:  

	mapgd ep -i metapopulation.pro -p 1 -o population1_allelfrequencies.txt
	mapgd cp -i metapopulation.pro -p 1,2 -o allelfrequency_comparison.txt

Alternatively, if the .pro file contains individual data you will want to run the allele command:

	mapgd allele -i metapopulation.pro -p 1 -o population1_allelfrequencies.txt

<h5> The above work flow using  I/O redirection </h5>
	samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam 
	| mapgd proview | mapgd allele -p 1,2 > allelefrequency.txt

or equivalently: 

	samtools mpileup -q 25 -Q 25 -B population1.sort.bam population2.sort.bam 
	| mapgd allele -p 1,2 > allelefrequency.txt

In the case where the allele command has been used the txt files can be further analized to estimate the seven IBD coefficients. Currently The tool mapgdutils.py must be used on the reformat the results, folowed by use of the relatedness.py script:

	mapgdutils.py -l 10 -e 0.01 -c 50 -C 200 --pro population.pro --map population.mpa > population.gcf
	convert2bcf population.gcf population.bcf
	relatedness.py population.bcf

<h5> Other Useful Programs </h5>

This program is intended for use with ".pro" described in the previous section. A slightly modified version of the program sam2pro, written by Bernhard Haubold, is included in this package and can be run by typing "mapgd proview". Programs for converting other file formats to ".pro" files are available in stand alone form http://guanine.evolbio.mpg.de/mlRho/

The workflow described above requires the programs bwa and samtools.

To download bwa please visit http://bio-bwa.sourceforge.net/

To download samtools please visit http://www.htslib.org/

<h5> For windows users </h5>

The windows binaries for V 2.1 are currently unavailable. You may try to compile the source code yourself, as we have tried to refrain from using any platform specific libraries, or you can send me (Matthew) an e-mail telling me to get the binaries up ASAP. 

If the binaries did exist, you could obtain them by clicking the "Download ZIP" button. You would then be prompted to save or open the file MAPGD-master.zip. Extract this file to the directory of your choice, which should create the new directory "MAPGD-master". Open this directory in windows explorer and click on the file "RunMeWin-32.bat" if you have a 32 bit version of windows or the "RunMeWin-64.bat" if you have a 64 bit versions of windows. If you don't know what version of windows you have, then just try clicking on both. You will be prompted to enter the name of the file you wish to analyze. You can type "test\test.pileup" to analyzea  test data set. To analyze your own data, just drag and drop the file into MAPGD-master folder, click on the RunMeWin-32.bat or RunMeWin-64.bat and type the filename of the file. The output will be saved to the file output.txt.   

mapgd can also be run from the command prompt, which can be accessed pressing the Windows logo key and r key simultaneously then typing "cmd" into menu which appears.

Although mapgd can be run by windows users, bwa, samtools, or other bioinformatic programs may not be available for window users.

<h5> Linux or Mac users </h5>

After clicking the "Download ZIP" button you will be prompted to save or open the file MAPGD-master.zip. Save this this file to the directory of your choice, then go to this director in a terminal, for example /home/Foo/Downloads/

Then type:

	unzip MAPGD.zip
	cd MAPGD-master/src
	make

The program can be installed for all users of a computer by typing:

	sudo make install

Scripts for both Linux and Mac users are present in the top level directory, or the program can be run by typing "mapgd ep -i FILENAME" where FILENAME is the a .pro file.

Some Apple users may not have developmental tools installed. To install them please type

	xcode --install

<h3> Bigred2 </h3>

*make sure to specify the number of threads to use with the ppn option* 
i.e., #PBS -l ppn=16

A note to IU users: Bigred2 has several different programming environments. While we hope to support use int he cray development environment in the future, currently support is only available in the gnu programming environment. To compile the code you will have to type "module rm PrgEnv-cray" then "module load PrgEnv-gnu". Additionally, usage of the relatedness python scripts requires that you type "module load boost"

<h3> Karst </h3>


module load boost
module load intel

<h3> References </h3>

Please cite the following paper (once it is written) when publishing results derived from this program:

Maruki, Takahiro, and Michael Lynch. "Genotype-Frequency Estimation from High-Throughput Sequencing Data." Genetics 201.2 (2015): 473-486.

For understanding the output of the relatedness.py command, it may be useful to read 

Ackerman, M. S., P. Johri, K. Spitze and M. Lynch, 2015  A general statistical model for coefficients of relatedness and its application to the analysis of population-genomic data. In prep. 

<h4> Copyright </h4>

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

For a copy of the GNU General Public License write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
