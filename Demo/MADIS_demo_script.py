# -*- coding: utf-8 -*-
"""
Created on Wed May 17 08:25:11 2023

@author: Jana
"""

# Demo for the MADis 
## Select nonancestral line from ref

import pandas as pd
import numpy as np
import MADis_230905 as md


# Input, user-defined
file_gen = open("L1_demo-Glyma19G120400.vcf","rt") # vcf file for genotype data input
gene_nam = "Glyma.19G120400"   # The name of the analyzed gene from the vcf file
res_nam = "L1_demo_res-Glyma19G120400"  # A name of the result file name
keep = 100  # The number of the best-scored combinations to save into the result file
max_com = 7 # The maximal number of the position in the analyzed combination  
change_line=37806160
file_phen = open("Pod_col_Soy1066.txt","rt")  # txt file with phenotype data


# Create a genotype matrix, change Alt 1 as ref for the selected line (in case of no MUT Ref sample with presumed ancestral Alt 1 allele at the pos)
geno=md.import_vcf_gene(file_gen,gene_nam,change_line=change_line)

# Create a phenotype vector and preselect the valid samples
nam=list()
phen=list()
for line in file_phen:
    line=line.split()
    if line[1]=="Bl":     # Define the WT phenotype
        nam.append(line[0])
        phen.append(0)
    elif line[1]=="Tn" or line[1]=="Br":     # Define the MUT phenotype
        nam.append(line[0])
        phen.append(1)
file_phen.close() 
print(f"Phen count {len(phen)}.")
wt_p=np.count_nonzero(np.array(phen)==0)
mut_p=np.count_nonzero(np.array(phen)==1) 
print(f"Phen count wt {wt_p}.")
print(f"Phen count mut {mut_p}.") 

# Adjust the genotype matrix to the phenotype vector (keep only the valid samples)
geno=geno[nam]
print(f"Genotype pos count {len(geno)}.")

# Run the MADIS analysis
res=pd.DataFrame(columns=["Pos_comb","N_pos_comb","Score","Score_max","N_WT_pheno","N_WT_corr_pred","N_MUT_pheno","N_MUT_corr_pred"])
for N in range(2,max_com+1):
    df=md.make_score(geno, phen, N)  # Run MADIS analysis for a selected number of positions in combination
    res = pd.concat([res,df])    # Combine the result for the different position combination count


# Save the analysis result
res=res.sort_values(by=["Score","N_pos_comb"], ascending=[False,True])
res=res.head(keep)
res.insert(0, 'Order', [a for a in range(1,len(res)+1)])
res.to_csv(f"{res_nam}.csv",index=False)     
