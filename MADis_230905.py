# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 09:58:40 2023

@author: biovja00
"""

import pandas as pd
import numpy as np
from itertools import combinations 


###############################################################################
#MADIS predefined functions 
###############################################################################


###############################################################################

# Transform vcf data into the binary matrix (part1.1), transform vcf line
def transform_line(line):
    ln=[int(line[1])]
    for a in line[9:]:
        if a[0]==a[2]:
            try:
              int(a[0])
              if a[0]=="0":ln.append(0)
              else:ln.append(1)
            except:ln.append(10)
        else: 
            ln.append(10)
            
                    
    return ln

# Transform vcf data into the binary matrix (part1.2), transform vcf line, and replace ref to alt
def transform_line_change(line):
    ln=[int(line[1])]
    for a in line[9:]:
             if a[0]==a[2]:
                 try:
                   int(a[0])
                   if a[0]=="1":ln.append(0)
                   else:ln.append(1)
                 except:ln.append(10)
             else: 
                 ln.append(10)
    return ln


#Transform vcf data into binary matrix (part 2.1), select lines for transformation, and collect transformed lines
### change_line: position of the line (int) for alt/ref change (None or position number for change (int))
### gene_nam: name of selected gene(must be included in the vcf file INFO column)
### eff: allele change predicted effect on protein (must be included in vcf file INFO column), (MOD or None)
def import_vcf_gene(file,gene_nam,eff="MOD",change_line=None):
    lst=list()
    n="n"
    for line in file:
        if n=="y":
            if gene_nam in line:
                line=line.split()
                if eff=="MOD":
                    if "HIGH" in line[7] or "MODERATE" in line[7]:  
                        if change_line==None:
                            line=transform_line(line)
                            lst.append(line)
                        else:
                            if line[1]==str(change_line):
                                line=transform_line_change(line)
                                lst.append(line)
                            else:
                                line=transform_line(line)
                                lst.append(line)
                else:
                    if change_line==None:
                            line=transform_line(line)
                            lst.append(line)
                    else:
                            if line[1]==str(change_line):
                                line=transform_line_change(line)
                                lst.append(line)
                            else:
                                line=transform_line(line)
                                lst.append(line)
        elif line.startswith("#CHROM"):
            hd=["POS"]+line.split()[9:]
            n="y"
    geno=pd.DataFrame(lst,columns=hd)
    geno=geno.set_index("POS")
    #print(df)
    file.close()
    return geno 

#Transform vcf data into binary matrix (part 2.2), select lines for transformation, and collect transformed lines
### change_line: position of the line (int) for alt/ref change (None or position number for change (int))
### chrom: chromosome value for the selected area (int)
### start|stop: position limit values for the selected area (int)
def import_vcf_pos(file,chrom,start,stop,change_line=None):
    lst=list()
    n="n"
    for line in file:
        if n=="y":
            line=line.split()
            ch=""
            for a in line[0]:
                if a.isdigit(): ch=ch+a
            if chrom==int(ch):
                if int(line[1])>=start and int(line[1])<=stop:

                    if change_line==None:
                            line=transform_line(line)
                            lst.append(line)
                    else:
                            if line[1]==str(change_line):
                                line=transform_line_change(line)
                                lst.append(line)
                            else:
                                line=transform_line(line)
                                lst.append(line)
                
        elif line.startswith("#CHROM"):
            hd=["POS"]+line.split()[9:]
            n="y"
    geno=pd.DataFrame(lst,columns=hd)
    geno=geno.set_index("POS")
    #print(df)
    file.close()
    return geno 

# MADis score calculation
d_score ={"10":-1,"01":-1,"20":-3,"21":-3,"30":-6,"31":-6,"40":-10,"41":-10,"50":-15,"51":-15,
            "60":-21,"61":-21,"70":-28,"71":-28,"80":-36,"81":-36,"90":-45,"91":-45}
def make_score(geno,phen,N,d_score=d_score):
    max_score=len(phen)
    wt_p=np.count_nonzero(np.array(phen)==0)
    mut_p=np.count_nonzero(np.array(phen)==1) 
    res=list()
    for com in combinations(geno.index, N):   #cout all combination of N for the position
        com=list(com)
        
        array=geno.loc[com]
        geno_sum=np.sum(array,axis = 0)
            
        geno_red=geno_sum%10
        
        score1,wt,mut=0,0,0
        for g,p in zip(geno_red,phen):
           
            if g==p: 
                score1=score1+1
                if p==0:wt=wt+1
                elif p==1:mut=mut+1
            else:
                score1=score1+d_score[str(g)+str(p)]
        N=N
        res.append([com,N,score1,max_score,wt_p,wt,mut_p,mut])  
            
    df=pd.DataFrame(res,columns=["Pos_comb","N_pos_comb","Score","Score_max","N_WT_pheno","N_WT_corr_pred","N_MUT_pheno","N_MUT_corr_pred"])
    df=df.sort_values(by=["Score"], ascending=False)  
    df["Explained_phen(%)"]=round(100/df["Score_max"]*(df["N_MUT_corr_pred"]+df["N_WT_corr_pred"]),2) 
   
    return df

###############################################################################