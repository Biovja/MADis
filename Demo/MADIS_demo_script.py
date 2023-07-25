# -*- coding: utf-8 -*-
"""
Created on Wed May 17 08:25:11 2023

@author: Jana
"""

import pandas as pd

import MADis_230725 as md


# Input
file_gen = open("Glyma_09G235100_gene_from_Soy1066.vcf","rt")
gene_nam = "Glyma.09G235100"
res_nam = "Glyma_09G235100_MADIS_res"
keep = 100
file_phen = open("Hilum_col_BlXBr_Soy1066_pheno_template.txt","rt")


# create a genotype matrix
geno=md.import_vcf(file_gen,gene_nam)

# create a phenotype vector
nam=list()
phen=list()
for line in file_phen:
    line=line.split()
    if line[1]=="0":
        nam.append(line[0])
        phen.append(0)
    elif line[1]=="1":
        nam.append(line[0])
        phen.append(1)
file_phen.close() 

# adjust the genotype matrix to the phenotype vector
geno=geno[nam]

# run MADIS analysis
res=pd.DataFrame(columns=["Pos_comb","N_pos_comb","Score","Score_max","N_WT_pheno","N_WT_corr_pred","N_MUT_pheno","N_MUT_corr_pred"])
for N in [2,3,4,5,6,7]:
    df=md.make_score(geno, phen, N)
    res = pd.concat([res,df])


# save analysis result
res=res.sort_values(by=["Score","N_pos_comb"], ascending=[False,True])
res=res.head(keep)
res.insert(0, 'Order', [a for a in range(1,len(res)+1)])
res.to_csv(f"{res_nam}.csv",index=False)     
