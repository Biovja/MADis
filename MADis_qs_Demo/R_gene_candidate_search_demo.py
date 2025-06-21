# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 08:28:02 2025

@author: Jana
"""

import pandas as pd
import MADis_qs as md
import numpy as np

# A name of the result file 
res_nam = "Demo_R_gene_scan_for_candidate_gene_in_loci.csv" 

# Input genotype, user-defined
file = open("Demo-candidate_search_genotype_R-gene.vcf", "rt") # vcf file for genotype data input 
### Create a genotype matrix
geno=md.import_vcf_pos(file,9,45765629-500000,45765629+500000)
file.close()
#print(geno)

# Input phenotype, user-defined
file_phen = open("Hilum_col_BlXBr_Soy1066.txt", "rt") # vcf file for genotype data input 
### Create a phenotype vector and preselect the samples with valid 
nam=list()
wtnam=list()
phen=list()
for line in file_phen:
    line=line.split()
    if line[1]=="Bl":     # Define the WT phenotype
        nam.append(line[0])
        wtnam.append(line[0])
        phen.append(0)
    elif line[1]=="Br":     # Define the MUT phenotype
        nam.append(line[0])
        phen.append(1)
file_phen.close() 
print(f"Phen count {len(phen)}.")
wt_p=np.count_nonzero(np.array(phen)==0)
mut_p=np.count_nonzero(np.array(phen)==1) 
print(f"Phen count wt {wt_p}.")
print(f"Phen count mut {mut_p}.")

# Data filtration 
### Adjust the genotype matrix to the phenotype vector (keep only the valid samples)
geno=geno[nam]
### 1st filtration. Remove the only ref pos for all valid samples
geno=geno[geno.sum(axis=1) > 0]
### 2nd filtration. Remove if too many samples with WT phen and Alt genotype
genowt=geno[wtnam]
genowt["Sum"]=genowt.sum(axis=1, numeric_only=True)
genowt=genowt[genowt["Sum"]<(len(wtnam)*0.1)]
keeppos=genowt.index.to_list()
geno=geno.loc[keeppos]

# Input gene info
genes_annot_tab=md.sel_genes("Chr09",45765629-500000,45765629+500000)
#print(genes_annot_tab)

# Analyze the genes for candidates
res=md.quick_candidates(geno,phen,genes_annot_tab)
print(res)

# Save the results
res.to_csv(f"{res_nam}.csv",index=False)
