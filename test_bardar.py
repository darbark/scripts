#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import os
import tarfile


# In[ ]:


# create dataframe from file in archive
def get_df(path, fname):
        with tarfile.open(f'{path}/{archive_name}') as t:
            with t.extractfile(f'Sample.kallisto_all/{fname}') as c:
                content = c.readlines()
                for si in range(len(content)):
                    content[si] = content[si].decode().rstrip('\n').split('\t')
                df = pd.DataFrame(content[1:], columns=content[0])
                return df
            
# compares three protocols, returns smallest / largest value            
def select_genes_with_lower_expr(selected, other1, other2, df):
    index_list = df[(df[selected] < df[other1]) & (df[selected] < df[other2])].index.tolist()
    low = set(index_list)
    return low

def select_genes_with_over_expr(selected, other1, other2, df):
    index_list = df[(df[selected] > df[other1]) & (df[selected] > df[other2])].index.tolist()
    high = set(index_list)
    return high

# find pathes to directories from df_annotation row
def find_path(df_annotation, path_to_start_dir):
    covaris = ''
    purigen = ''
    qiagen = ''
    for root, dirs, _ in os.walk(path_to_start_dir):
        if df_annotation.iloc[i, 2] in dirs:
            covaris = os.path.join(root, df_annotation.iloc[i, 2], 'output')
        if df_annotation.iloc[i, 3] in dirs:
            purigen = os.path.join(root, df_annotation.iloc[i, 3], 'output')  
        if df_annotation.iloc[i, 4] in dirs:
            qiagen = os.path.join(root, df_annotation.iloc[i, 4], 'output')
        if covaris and purigen and qiagen:
            break
    return covaris, purigen, qiagen

def upgrade_df(df):
    for col in df.columns.values.tolist()[1:]:
        df[col] = pd.to_numeric(df[col])
    df['log_norm'] = np.log2((df['tpm'] + 1))
    return df

def columns_for_gene_expression_compare(low, high, codegenes):
    low_num = len(low[annotation_trunc.iloc[i, 0]])
    low_cod = round((low_num / len(codegenes)), 2)
    high_num = len(high[annotation_trunc.iloc[i, 0]])
    high_cod = round((high_num / len(codegenes)), 2)
    data_list = [annotation_trunc.iloc[i, 0], annotation_trunc.iloc[i, 1], low_num, low_cod, high_num, high_cod]
    return data_list
            


# In[ ]:


# variables with used files
archive_name = 'expression-kallisto-xena.tar.gz'
logfile = 'Sample-kallisto-Xena-gene-log2TPMplus1.tsv'
initialfile = 'abundance.tsv'
withoutnc = 'Sample-kallisto-Xena-gene-log2TPMplus1_without_noncoding.tsv'

# absolute path to the directory 'test_task'
outside_path = '/home/dasha/Documents/test_task'


# In[ ]:


annotation = pd.read_csv(f'{outside_path}/FFPE_annotation.txt', sep='\t')
annotation


# In[ ]:


sns.set_style("darkgrid")
for i in range(len(annotation)):
    covaris_path, purigen_path, qiagen_path = find_path(annotation, outside_path)
        
    # create the dataframes from log2TPMplus1 file and use it to plot kde graph 
    
    fig, (ax1, ax2, ax3) = plt.subplots(figsize=(18, 5), ncols=3, sharey=True)
    plt.suptitle(f'{annotation.iloc[i, 1]} (sample {annotation.iloc[i, 0]})', fontsize=22)
    
    if covaris_path:
        covaris_df = get_df(covaris_path, logfile)
        sns.kdeplot(data=covaris_df['Sample'], label='covaris', shade=True, color='orange', ax=ax1)
        ax1.legend(fontsize='xx-large')
        ax1.set(xlim=(0, None))
                     
    if purigen_path:
        purigen_df = get_df(purigen_path, logfile)
        sns.kdeplot(data=purigen_df['Sample'], label='purigen', shade=True, color='lawngreen', ax=ax2)
        ax2.legend(fontsize='xx-large')
        ax2.set(xlim=(0, None))
                      
    if qiagen_path:
        qiagen_df = get_df(qiagen_path, logfile)
        sns.kdeplot(data=qiagen_df['Sample'], label='qiagen', shade=True, color='darkturquoise', ax=ax3)
        ax3.legend(fontsize='xx-large')
        ax3.set(xlim=(0, None))


# In[ ]:


c_low = dict()
c_high = dict()

p_low = dict()
p_high = dict()

q_low = dict()
q_high = dict()

# use only rows with all data
annotation_trunc = annotation.replace('Failed', np.nan).dropna()

for_covaris = list()
for_purigen = list()
for_qiagen = list()

for i in range(len(annotation_trunc)):
    covaris_path, purigen_path, qiagen_path = find_path(annotation_trunc, outside_path)
            
    if covaris_path:
        covaris_df = get_df(covaris_path, logfile).set_index('Gene').rename(columns={'Sample': 'Covaris'})
        covaris_coding_genes = set(get_df(covaris_path, withoutnc)['Gene'].tolist())
        
    if purigen_path:
        purigen_df = get_df(purigen_path, logfile).set_index('Gene').rename(columns={'Sample': 'Purigen'})
        purigen_coding_genes = set(get_df(purigen_path, withoutnc)['Gene'].tolist())
               
    if qiagen_path:
        qiagen_df = get_df(qiagen_path, logfile).set_index('Gene').rename(columns={'Sample': 'Qiagen'})
        qiagen_coding_genes = set(get_df(qiagen_path, withoutnc)['Gene'].tolist())
    
    # merge data by genes names and keep only rows with protocol-differ expression     
    compare = pd.concat([covaris_df, purigen_df, qiagen_df], axis=1)
    for col in compare.columns.values.tolist():
        compare[col] = pd.to_numeric(compare[col])
    compare = compare.round(0)
    inames = compare[(compare['Covaris'] == compare['Purigen']) & (compare['Purigen'] == compare['Qiagen'])].index
    compare.drop(inames , inplace=True)    
    
    # for each sample collect a set of genes with differ expression
    c_low[annotation_trunc.iloc[i, 0]] = select_genes_with_lower_expr('Covaris', 'Purigen', 'Qiagen', compare)
    c_high[annotation_trunc.iloc[i, 0]] = select_genes_with_over_expr('Covaris', 'Purigen', 'Qiagen', compare)
    
    p_low[annotation_trunc.iloc[i, 0]] = select_genes_with_lower_expr('Purigen', 'Covaris', 'Qiagen', compare)
    p_high[annotation_trunc.iloc[i, 0]] = select_genes_with_over_expr('Purigen', 'Covaris', 'Qiagen', compare)
    
    q_low[annotation_trunc.iloc[i, 0]] = select_genes_with_lower_expr('Qiagen', 'Purigen', 'Covaris', compare)
    q_high[annotation_trunc.iloc[i, 0]] = select_genes_with_over_expr('Qiagen', 'Purigen', 'Covaris', compare)
   
    # save rows for dataframes
    for_covaris.append(columns_for_gene_expression_compare(c_low, c_high, covaris_coding_genes))
    for_purigen.append(columns_for_gene_expression_compare(p_low, p_high, purigen_coding_genes))
    for_qiagen.append(columns_for_gene_expression_compare(q_low, q_high, qiagen_coding_genes))
    
columns = ['Sample', 'Tissue', 'Low_num', 'Low_coding', 'Over_num', 'Over_coding']

covaris_shift = pd.DataFrame(for_covaris, columns=columns)
purigen_shift = pd.DataFrame(for_purigen, columns=columns)
qiagen_shift = pd.DataFrame(for_qiagen, columns=columns)


# In[ ]:


# ##### Контент колонок: 
# * **Low_num** - количество генов, экспрессия которых при данной пробоподготовке ниже, чем при других
# * **Low_coding** - доля кодирующих генов в *Low_num*
# * **Over_num** - количество генов, экспрессия которых при данной пробоподготовке выше, чем при других
# * **Over_coding** - доля кодирующих генов в *Over_num*


# In[ ]:





# In[ ]:


covaris_shift


# In[ ]:


purigen_shift


# In[ ]:


qiagen_shift


# In[ ]:


sns.set_style("darkgrid")
x = 'length'
y = 'log_norm'
# ylim = (0.4, 15) 
# xlim = (5000,100000)

# create the dataframes from abundance file, calculate additional column with log(tpm + 1) normalisation
# and use it to plot graph 

for i in range(len(annotation)):
    covaris_path, purigen_path, qiagen_path = find_path(annotation, outside_path)
        
    if covaris_path:
        covaris_df = get_df(covaris_path, initialfile)
        covaris_df = upgrade_df(covaris_df)
        
        sns.jointplot(x=x, y=y, data=covaris_df, kind='scatter', color='orange')
        plt.title(f'Covaris, {annotation.iloc[i, 1]} (sample {annotation.iloc[i, 0]})', y=1.3, fontsize = 16)

    
    if purigen_path:
        purigen_df = get_df(purigen_path, initialfile)
        purigen_df = upgrade_df(purigen_df)
        
        sns.jointplot(x=x, y=y, data=purigen_df, kind='scatter', color='lawngreen')
        plt.title(f'Purigen, {annotation.iloc[i, 1]} (sample {annotation.iloc[i, 0]})', y=1.3, fontsize = 16)

        
    if qiagen_path:
        qiagen_df = get_df(qiagen_path, initialfile)
        qiagen_df = upgrade_df(qiagen_df)
        
        sns.jointplot(x=x, y=y, data=qiagen_df, kind='scatter', color='darkturquoise')
        plt.title(f'Quagen, {annotation.iloc[i, 1]} (sample {annotation.iloc[i, 0]})', y=1.3, fontsize = 16)

