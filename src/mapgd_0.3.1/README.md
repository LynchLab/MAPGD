<script src='https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'></script>

<h1> MAPGD version 0.3.1 </h1>

<h3>Contents </h3>

####[Introduction](https://lynchlab.github.io/MAPGD/index.html#-introduction-)
####[Basic Design](https://lynchlab.github.io/MAPGD/index.html#-basic-design-)
####[Tutorials](https://lynchlab.github.io/MAPGD/index.html#-tutorials-)

<h3> Introduction </h3>

MAPGD is a series of related programs that estimate allele frequency, heterozygosity, Hardy-Weinberg disequilibrium and identity by descent (IBD) coefficients from population genomic data using statistically rigorous maximum likelihood approach.   

We have tried to produce an extensible program ...

Ultimately the genetic structure of a population is fully specified by the genotypes of the individuals that compose that popualtion. This means that if we can accurately calculate all [genotypic probabilities](https://lynchlab.github.io/MAPGD/classgenotype.html), then the calculation of any other population statistics becomes trivial. However, in order to calcule genotypic probabilities we must take acount of the errors made in the genotyping process, include infering the precense of alleles that are not there, which can arrise from sequencing error or mistakes made aligning to a reference, and failing to detect the precense of alleles which are there because of low coveraged or biased sequencing of a single parental chromosome. We maximize a likelihood equation to account for sequeincing error and failure to sample genotypes, and then we test the data for fit to the parameters, and reject estimates where the data has a poor fit to the estimated parameters. 

<h3> Bais Design </h3>

$ P( \boldsymbol \theta | \boldsymbol X  )= \frac {P( \boldsymbol X | \boldsymbol \theta  ) P( \boldsymbol \theta  ) } { P( \boldsymbol X  ) } $	

<h4> Style Guidlines </h4>
mapgd is written to conform to the [GNU style guidlines](www.gnu.org/prep/standards/standards.html#Formatting), at least to the extent that I have had time to read and implement the guidlines. 
Class and type names should be whole english words starting with an initial capital (to distinguish them from the core classes and types of c++). We do not use CamelCase. Variable names should be all lower case with words seperated by an underscore. Private variables ...

* Internationalization: These are good standards... Not implemented yet. 

<h4> Likelihood models:</h4>

<h5> Pior:</h5>

<h5> Likelihoods:</h5>

<h5> Posteriors:</h5>

###Maximizing Priors:

###Genotypic Posteriors:

###Summary Statistics:

####[MAP files](https://lynchlab.github.io/MAPGD/classmap__file.html)
####[PRO files](https://lynchlab.github.io/MAPGD/classmap__file.html)

<h3> Tutorials </h3>

* [An introduction to quartets](https://lynchlab.github.io/MAPGD/tutorial/quartet.md)
* [Reading and writing to map files](https://lynchlab.github.io/MAPGD/tutorial/map.md)
* [Making likelihood functions](https://lynchlab.github.io/MAPGD/tutorial/likelihood.md)
* [Maximizing a likelihood function](https://lynchlab.github.io/MAPGD/tutorial/maximize.md)
