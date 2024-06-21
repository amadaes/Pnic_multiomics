### Python3 script that fixes the MSTRG.x ids added by stringtie and combines them with a count matrix.
### It uses the output from stringtie merge - file2 and a count matrix to generate a single file with the reference IDs and counts
### An absolute path to the files can be specified in the path variable or individually for file1 and file2
### With some changes, the script can also be used to only separate and select specific features from a gtf/ gff file
### Run script with Python3 merge_stringtie_featurecounts.py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os

#input path 
path='D:/Activities/DESeq/RAW_DATA/all_samples'
os.chdir(path)

file1=path+'/stringtie2/featureCounts/counts_gene.txt'
df=pd.read_csv(file1,sep='\t',index_col=0, header=1)
df.index.rename('gene_id', inplace=True) #rename index

file2=path+'/stringtie/stringtie.merged.gtf'
df2=pd.read_csv(file2, sep='\t', comment='#',header=None)

features_column=8
def add_feature(df2,feature):
    df2[feature]=''
    for i in df2.index:
        g= df2.loc[i,features_column]
        try: 
            temp=np.array(g.split('; '))
            temp2=np.array([t.split(' ')[1] if feature == t.split(' ')[0] else ''  for t in temp])
            temp3=temp2[temp2!=''][0].replace("\"",'')
            temp3=temp3.replace(" ","").replace(";","")
            df2.loc[i,feature]=temp3
        except:
            1
            #print ("missing at position "+ str(i))

#these functions add one by one the features of interest; can be changed with whatever is of interest
add_feature(df2,feature='ref_gene_id')
add_feature(df2,feature='gene_id')
add_feature(df2,feature='gene_name')

#set index of df2 to stringtie ids so they match index of featureCounts df
df2.set_index('gene_id',drop=0,inplace=True)

#this selects only transcript rows; can be changed to exon, gene etc.  
df2=df2[df2[2]=="transcript"][[0,3,4,"ref_gene_id", "gene_name"]]

#rename start and end columns
df2=df2.rename({0:'Chr',3:'Start',4:'End'}, axis=1) 

#combine tables into a new dataframe raw_counts_ids using the gene_id, Start and End coordinates
raw_counts_ids=pd.merge(df,df2,how='left', on=['gene_id','Start','End'], indicator=True)

#fill in the ref_gene_id which are not available from ref database with MSTRG IDs from stringtie merge (index keys)
raw_counts_ids['ref_gene_id'] = raw_counts_ids['ref_gene_id'].replace('', pd.NA).fillna(raw_counts_ids.index.to_series())

#fill in the gene_name which are not available from ref database with ref_gene_id (now including MSTRG IDs from stringtie merge)
raw_counts_ids['gene_name'] = raw_counts_ids['gene_name'].replace('', pd.NA).fillna(raw_counts_ids['ref_gene_id'])

#make ref_gene_id first column
first_column=raw_counts_ids.pop('ref_gene_id')
raw_counts_ids.insert(0,'ref_gene_id', first_column)

#make gene_name the second column
second_column=raw_counts_ids.pop('gene_name')
raw_counts_ids.insert(1,'gene_name', second_column)
df3=pd.DataFrame(raw_counts_ids)

#clean up the table to exclude redundant info or unnecessary for next step
#df3=df3.drop('gene_name', axis=1)
df3=df3.drop('Chr_y', axis=1)
df3=df3.drop('_merge', axis=1)

#rename columns to match the file that the R script expects 
df3=df3.rename({'ref_gene_id':'Geneid','gene_name':'Genename','Chr_x':'Chr'}, axis=1) 

#save final data to csv file
print(df3.to_csv(path+'/ids_and_counts.txt', index=False))
