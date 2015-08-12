##MAPGD version 0.3.1

[<b>Download MAPGD.zip</b>](https://github.com/LynchLab/MAPGD/archive/master.zip) [(C) Michael Lynch](https://github.com/LynchLab/MAPGD#-copyright-)

[<b>Visit the MAPGD development page</b>](https://lynchlab.github.io/MAPGD/)

###Contents 
####[Introduction](https://github.com/LynchLab/MAPGD#-introduction-)
#####[FAQ](https://github.com/LynchLab/MAPGD#-faq-)
#####[Changes](https://github.com/LynchLab/MAPGD#-changes-)
####Basic Usage
#####[Commands](https://github.com/LynchLab/MAPGD#-commands-)
#####[Input/Output](https://github.com/LynchLab/MAPGD#-inputoutput-)
#####[Scripting](https://github.com/LynchLab/MAPGD#-scripting-)
#####[Examples](https://github.com/LynchLab/MAPGD#-examples-)
#####[Linux or Mac users](https://github.com/LynchLab/MAPGD#-linux-or-mac-users-)
####Misc.
#####[Extras](https://github.com/LynchLab/MAPGD#-extras-)
#####[Bigred2](https://github.com/LynchLab/MAPGD#-bigred2-)
#####[References](https://github.com/LynchLab/MAPGD#-references-)

<h3> Introduction </h3>

MAPGD is a series of related programs that estimate allele frequency, heterozygosity, Hardy-Weinberg disequilibrium and identity by descent (IBD) coefficients from population genomic data using statistically rigorous maximum likelihood approach. It is primarily useful for the analysis of low coverage population genomic data, and provides minimum MSE estimators of these statistics on low coverage sequence. Although other tools, such as vcftools, give similar allele frequency estimates with less computational investment when coverage is high, MAPGD should always be preferred to other programs when calculating IBD coefficients. 

<h5> FAQ </h5>

<b> What does mapgd do? </b>

* estimate frequencies
* compare frequencies
* estimate genotypes
* compare genotypes

<b> Why don't you use a vcf format? </b>

We intend to implement read/writing in the vcf/bcf format as quickly as possible to increase the compatibility of our program with existing tools. 

<b> Why don't you provide information on indels? </b>

We currently do not have a likelihood model to account for errors in calling indels, so we cannot incorporate indels into our program at this time. This is on our TODO list, but may not occur for some time.

<b> How long does it take to run? </b>

Typical benchmarks with 16 threads on a 2.6 GHz put us at around 18,000 sites a second for 96 simulated individuals at 10x coverage. This means that the typical invertebrate population will take around three hours to analyze on a good computer, and a vertebrate genome might take a few days. If you have managed to sequence *Paris japonica* you're looking at three months of computation time. We are not currently focused on scaling our program up for larger genomes, however, if you want to analyze something like *japonica*, let us know and we will make MPS and CUDA priorities! 

<b> Help, I can't get the program to compile. </b>

MAPGD has been written using the gnu C++11 style, but I've attempted to minimize my reliance on the newer features of the C++ in the hopes that older styles of C++ will have greater cross-platfrom compatibility. However, I have not yet made a systematic examination of mapgd's cross-platform copmpatibility. I have compiled mapgd on :

* Ubuntu Linux, 14.04

* Red Hat Enterprise Linux 6.x

* A version of the Cray Linux Environment.

* OS X Yosemite.

To compile on OS X, you need to type 'make noomp' because the default OS X does compiler does not support openmp. I sould have this fixed shortly, after which OS X users can get back to typing 'make' like eveyone else, but it may be a few weeks before that makes it to the top of the TODO list. 
 
Version 0.1 will compile on windows 32 and 64 bit versions IIRC, but getting current versions of mapgd working on windows is near the bottom of the TODO list.

Long story short: If you have problems, please please please e-mail me, and I will work to get it to compile for you as quickly as I can. 

<b> Help, the program keeps crashing/hanging </b>
 
The first thing you should do is e-mail me (for contact infromation type 'mapgd -h'). There are a few known issues that may not be properly fixed in the near future. First off, mapgd has problems merging ddifferent mpileup files into a single .pro file when some of the mpileup files have missing scaffolds. Ultimately we hope to address this by reading directly from bam files which have headers containing scaffold name and size informatoin, however, until this is implemented you can print the header with the 'samtools view -H FILENAME.bam > FILENAME-header.bam' command, and then give this header to mapgd with the 'mapgd proview -H FILENAME-header.bam -i \*.mpileup' command. 

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

	git commit -m "SOMETHING MEANINGFULL"

Finally show me your changes by typing.

	git push

<h5> Changes </h5>

Data is now saved to a sqlite3 database by default. This data can be viewed using a wide variety of software, e.g [sqlitebrowser](http://sqlitebrowser.org/). This format is our stable format we guarand useful than formats used in other software. Also, we have had an extensive overhaul of the basic usage, so please re-read that section.

![SQL screenshots](https://github.com/LynchLab/mapgd/docs/sql.png) 

<h5> Basic Examples </h5>

<b> Estimating Allele Frequency </b>:
<b> Cleaning the data </b>:
<b> Estimating Heterozygousity </b>:
<b> Estimating Silent Site Heterozygousity </b>:
<b> using SQL commands </b>

<h5> Commands </h5>

<b> read </b>:
<b> write </b>:
<b> estimate </b>:
<b> compare </b>:
<b> filter </b>:

<h4> Extras </h4>

<h6> mapgdutils.py </h6>The script mapgdutils.py is a stop gap measure used to generate genotypic likelihoods (i.e. gcf files) for the analysis of the relatedness of a pair of individuals. A script of this name is likely to persist in future version of the program to perform computationally simple task that may be useful; however, we hope to generation genotypic likelihoods directly in mapgd in the near future.

	mapgd proview -i \*.mpileups | mapgd ei -p FILENAME.pro -o FILENAME.map 
	python mapgdutils.py FILENAME.pro FILENAME.map > GENOTYPES.gcf

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

Currently mapgd allows for the estimation of allele frequencies in individually labeled and pooled population sequence through the “ei” and "ep" command respectively, the comparison of allele frequencies between two populations with the "cp" command,  the conversion of mpileups to the ".pro" format though the proview command, and the estimatation of seven components of pair-wise relatedness through the "rel" command.

For example, if the sequencing center gives you two files called "seq1.fastq" and "seq2.fastq" your entire work flow might look something like this: 

<h5> Other Useful Programs </h5>

This program is intended for use with ".pro" described in the previous section. A slightly modified version of the program sam2pro, written by Bernhard Haubold, is included in this package and can be run by typing "mapgd proview". Programs for converting other file formats to ".pro" files are available in stand alone form http://guanine.evolbio.mpg.de/mlRho/

The workflow described above requires the programs bwa and samtools.

To download bwa please visit http://bio-bwa.sourceforge.net/

To download samtools please visit http://www.htslib.org/

<h5> Linux or Mac users </h5>

After clicking the "Download ZIP" button you will be prompted to save or open the file MAPGD-master.zip. Save this this file to the directory of your choice, then go to this director in a terminal, for example /home/Foo/Downloads/

Then type:

	unzip MAPGD.zip
	cd MAPGD-master/src
	make

The program can be installed for all users of a computer by typing:

	sudo make install

Scripts for both Linux and Mac users are present in the top level directory, or the program can be run by typing "mapgd ep -i FILENAME" where FILENAME is the a .pro file.

<h3> Bigred2 </h3>

A note to IU users: Bigred2 has several different programming environments. While we hope to support use int he cray development environment in the future, currently support is only available in the gnu programming environment. To compile the code you will have to type "module rm PrgEnv-cray" then "module load PrgEnv-gnu".

<h3> References </h3>

Please cite the following paper (once it is written) when publishing results derived from this program:

? A Paper ?. 

For understanding the output of the relatedness.py command, it may be useful to read 

Ackerman, M. S., P. Johri, K. Spitze and M. Lynch, 2015  A general statistical model for coefficients of relatedness and its application to the analysis of population-genomic data. In prep. 

<h4> Copyright </h4>

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

For a copy of the GNU General Public License write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
