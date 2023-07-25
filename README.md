# MADis
MADIS is a Python-based tool designed to help with the identification of potential multiple allele CMs in the genes of interest. 
As multiple CMs are hard to discover with classical genomic tools like GWAS. The MADIS package uses user-provided phenotype and 
genotype data for the estimation of the potential multiple allele CMs. The correspondence between the set hypothesis and the 
tested combination of genotype positions and phenotype is expressed by scoring. The score determined the performance of the 
combination of genotype positions and phenotype as potential CMs for the identified trait. The MADIS results are the tab-based 
list of the best fitting genotype positions combination for the phenotype of interest with additional analysis information for 
easier evaluation of the provided results. The MADIS python package consists of 2 major functions for the user to utilize. The first
function import_vcf(file,gene_nam) enables to import of genotype data and transforms them into the genotype binary matrix used for MADIS
scoring. The second main function make_score(geno,phen,N,d_score=d_score) produces the result table with the potential multiple CMs scores 
and additional information.

The tool is developed in Python 3.8 programing language and utilizes the following Python packages: pandas(1.3.5, https://pandas.pydata.org/ accessible on 20.07.2023), 
NumPy(2.21.5, https://numpy.org/ accessed on 20.07.2023) and itertools(3.11.4, https://docs.python.org/3/library/itertools.html accessed on 20.07.2023).
