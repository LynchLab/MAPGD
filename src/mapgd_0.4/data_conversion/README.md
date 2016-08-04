Data conversion {#Conversion}
===========================

These classes provide the basic functionality to convert between data in 
external formats, such as sam, bed, or vcf, and data in internal formats.

In general external formats discard data critical for statistical analysis and
cannot be converted directly into files usable by mapgd. 

To perform these conversions we leverage the basic data-processing schema of 
mapgd, which requires that different 

Thus

GFF:
BED:
SAM→ PRO
PRO→ MAP
(MAP, PRO)→ GCF
GCF→ VCF
(MAP, PRO)→ LNK
GCF→ REL
