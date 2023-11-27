#!/usr/bin/env python

### Import packages ###
import numpy as np
import pandas as pd
import sys


### Import data frame for each sample and concatenate in one data frame ###
df_concat = pd.DataFrame()

for i in range(1, len(sys.argv)-3):
    df = pd.read_csv(sys.argv[i], sep="\t", na_values='.')
    n = sys.argv[i].find('.')
    name = sys.argv[i][0:n]
    df["barcode"] = name
    df_concat = pd.concat([df_concat, df], axis =0)


### Import file with the type of variants ###
variants_type = pd.read_csv(sys.argv[len(sys.argv)-3], sep="\t")

sample_names = pd.read_csv(sys.argv[len(sys.argv)-2], sep=",")

exons = pd.read_csv(sys.argv[len(sys.argv)-1], sep=",")


### Reset index ###
df_concat = df_concat.reset_index(drop = True)


### Rename columns ###
df_concat.columns = ['chromosome', 'position', 'ref', 'alt',
                     'allele_depth', 'sample_depth', 'variant_type', 'gene_id',
                     'transcript_id', 'HGVS.c', 'HGVS.p', 'exon_rank', 'non_cancer_AF',
                     'non_cancer_AF_eas', 'non_cancer_AF_asj', 'non_cancer_AF_afr',
                     'non_cancer_AF_amr', 'non_cancer_AF_nfe', 'non_cancer_AF_fin',
                     'non_cancer_AF_sas', 'cadd_raw_v1.6', 'cadd_phred_v1.6',
                     'cadd_raw_v1.4', 'cadd_phred_v1.4', 'VEST4_score', 'VEST4_rankscore', 
                     'MetaSVM_score', 'MetaSVM_rankscore','MetaSVM_pred', 'MetaLR_score', 
                     'MetaLR_rankscore', 'MetaLR_pred','Reliability_index_Meta_score', 
                     'REVEL_score', 'REVEL_rankscore','BayesDel_addAF_score', 
                     'BayesDel_addAF_rankscore','BayesDel_addAF_pred', 'BayesDel_noAF_score', 
                     'BayesDel_noAF_rankscore','BayesDel_noAF_pred', 'ClinPred_score', 
                     'ClinPred_rankscore','ClinPred_pred', 'DANN_score', 'DANN_rankscore', 
                     'clinvar_id','clinvar_clnsig', 'clinvar_trait', 'clinvar_review', 
                     'clinvar_hgvs','clinvar_var_source', 'clinvar_MedGen_id', 
                     'clinvar_OMIM_id','clinvar_Orphanet_id', 'barcode']


### Split depth for allele 1 and allele 2 ###
df_concat[['allele1_depth', 'allele2_depth']] = df_concat["allele_depth"].str.split(",",expand=True)


### Convert allele depth to numeric value ###
df_concat['allele1_depth'] = pd.to_numeric(df_concat['allele1_depth'])
df_concat['allele2_depth'] = pd.to_numeric(df_concat['allele2_depth'])


### Calculate the allelic ratio ###
df_concat['allelic_frequency'] = ''
for i in range(len(df_concat)):
    x = df_concat['allele1_depth'][i]
    y = df_concat['allele2_depth'][i]
    z = y / (x + y)
    df_concat['allelic_frequency'][i] = z.round(decimals = 2)


### Add ID information ###
df_concat = df_concat.merge(sample_names, how="inner", left_on="barcode", right_on="run_barcode")
df_concat = df_concat.drop("run_barcode", axis = 1)


### Define variables for filter variants  ###
AF_MIN = 0.2
AF_MAX = 0.8
DP = 10
MAF = 0.0005


### Filter by allele frequency (AF) and read depth (DP) ###
df_AF_DP = df_concat.loc[(df_concat["allelic_frequency"] >= AF_MIN)
                           & (df_concat["allelic_frequency"] <= AF_MAX)
                           & (df_concat["sample_depth"] >= DP)]


### Filter by minor allelic frequency (MAF) ###
df_MAF = df_AF_DP.loc[(df_AF_DP["non_cancer_AF"] <= MAF)
                      & (df_AF_DP["non_cancer_AF_eas"] <= MAF)
                      & (df_AF_DP["non_cancer_AF_asj"] <= MAF)
                      & (df_AF_DP["non_cancer_AF_afr"] <= MAF)
                      & (df_AF_DP["non_cancer_AF_amr"] <= MAF)
                      & (df_AF_DP["non_cancer_AF_nfe"] <= MAF)
                      & (df_AF_DP["non_cancer_AF_sas"] <= MAF)
                      | (df_AF_DP["non_cancer_AF"].isna())
                      | (df_AF_DP["non_cancer_AF_eas"].isna())
                      | (df_AF_DP["non_cancer_AF_asj"].isna())
                      | (df_AF_DP["non_cancer_AF_afr"].isna())
                      | (df_AF_DP["non_cancer_AF_amr"].isna())
                      | (df_AF_DP["non_cancer_AF_nfe"].isna())
                      | (df_AF_DP["non_cancer_AF_sas"].isna())]


### Round MAF ###
df_MAF['non_cancer_AF']     = df_MAF['non_cancer_AF'].round(decimals = 5)
df_MAF['non_cancer_AF_eas'] = df_MAF['non_cancer_AF_eas'].round(decimals = 5)
df_MAF['non_cancer_AF_asj'] = df_MAF['non_cancer_AF_asj'].round(decimals = 5)
df_MAF['non_cancer_AF_afr'] = df_MAF['non_cancer_AF_afr'].round(decimals = 5)
df_MAF['non_cancer_AF_amr'] = df_MAF['non_cancer_AF_amr'].round(decimals = 5)
df_MAF['non_cancer_AF_nfe'] = df_MAF['non_cancer_AF_nfe'].round(decimals = 5)
df_MAF['non_cancer_AF_fin'] = df_MAF['non_cancer_AF_fin'].round(decimals = 5)
df_MAF['non_cancer_AF_sas'] = df_MAF['non_cancer_AF_sas'].round(decimals = 5)


### Reset index ### 
df_MAF = df_MAF.reset_index(drop = True)


### Round CADD raw score ###
df_MAF['cadd_raw_v1.6'] = df_MAF['cadd_raw_v1.6'].round(decimals = 3)
df_MAF['cadd_raw_v1.4'] = df_MAF['cadd_raw_v1.4'].round(decimals = 3)


### Select the first VEST4 score ###
for i in range(len(df_MAF)):
    if isinstance(df_MAF['VEST4_score'][i], str):
        n = df_MAF['VEST4_score'][i].find(',')
        if n > 0:
            df_MAF['VEST4_score'][i] = df_MAF['VEST4_score'][i][0:n]
        elif n < 0:
            df_MAF['VEST4_score'][i] = df_MAF['VEST4_score'][i]


### Round VEST4 score ###
df_MAF['VEST4_score'] = pd.to_numeric(df_MAF['VEST4_score'], errors = 'coerce')
df_MAF['VEST4_score'] = df_MAF['VEST4_score'].round(decimals = 3)


### Select the first VEST4 rankscore ###
for i in range(len(df_MAF)):
    if isinstance(df_MAF['VEST4_rankscore'][i], str):
        n = df_MAF['VEST4_rankscore'][i].find(',')
        if n > 0:
            df_MAF['VEST4_rankscore'][i] = df_MAF['VEST4_rankscore'][i][0:n]
        elif n < 0:
            df_MAF['VEST4_rankscore'][i] = df_MAF['VEST4_rankscore'][i]


### Round VEST4 rankscore ###
df_MAF['VEST4_rankscore'] = pd.to_numeric(df_MAF['VEST4_rankscore'], errors = 'coerce')
df_MAF['VEST4_rankscore'] = df_MAF['VEST4_rankscore'].round(decimals = 3)


### Select first META SVM score ###
for i in range(len(df_MAF)):
    if isinstance(df_MAF['MetaSVM_score'][i], str):
        n = df_MAF['MetaSVM_score'][i].find(',')
        if n > 0:
            df_MAF['MetaSVM_score'][i] = df_MAF['MetaSVM_score'][i][0:n]
        elif n < 0:
            df_MAF['MetaSVM_score'][i] = df_MAF['MetaSVM_score'][i]


### Round META SVM score ###   
df_MAF['MetaSVM_score'] = pd.to_numeric(df_MAF['MetaSVM_score'], errors ='coerce')
df_MAF['MetaSVM_score'] = df_MAF['MetaSVM_score'].round(decimals = 3)


### Select first META SVM rankscore ###
for i in range(len(df_MAF)):
    if isinstance(df_MAF['MetaSVM_rankscore'][i], str):
        n = df_MAF['MetaSVM_rankscore'][i].find(',')
        if n > 0:
            df_MAF['MetaSVM_rankscore'][i] = df_MAF['MetaSVM_rankscore'][i][0:n]
        elif n < 0:
            df_MAF['MetaSVM_rankscore'][i] = df_MAF['MetaSVM_rankscore'][i]


### Round META SVM rankscore ###
df_MAF['MetaSVM_rankscore'] = pd.to_numeric(df_MAF['MetaSVM_rankscore'], errors = 'coerce')
df_MAF['MetaSVM_rankscore'] = df_MAF['MetaSVM_rankscore'].round(decimals = 3)


### Select first META SVM pred ###
for i in range(len(df_MAF)):
    if isinstance(df_MAF['MetaSVM_pred'][i], str):
        n = df_MAF['MetaSVM_pred'][i].find(',')
        if n > 0:
            df_MAF['MetaSVM_pred'][i] = df_MAF['MetaSVM_pred'][i][0:n]
        elif n < 0:
            df_MAF['MetaSVM_pred'][i] = df_MAF['MetaSVM_pred'][i]


### Select first META LR score ### 
for i in range(len(df_MAF)):
    if isinstance(df_MAF['MetaLR_score'][i], str):
        n = df_MAF['MetaLR_score'][i].find(',')
        if n > 0:
            df_MAF['MetaLR_score'][i] = df_MAF['MetaLR_score'][i][0:n]
        elif n < 0:
            df_MAF['MetaLR_score'][i] = df_MAF['MetaLR_score'][i]


### Round META LR score ###
df_MAF['MetaLR_score'] = pd.to_numeric(df_MAF['MetaLR_score'], errors = 'coerce')
df_MAF['MetaLR_score'] = df_MAF['MetaLR_score'].round(decimals = 3)


### Select first META LR rankscore ###
for i in range(len(df_MAF)):
    if isinstance(df_MAF['MetaLR_rankscore'][i], str):
        n = df_MAF['MetaLR_rankscore'][i].find(',')
        if n > 0:
            df_MAF['MetaLR_rankscore'][i] = df_MAF['MetaLR_rankscore'][i][0:n]
        elif n < 0:
            df_MAF['MetaLR_rankscore'][i] = df_MAF['MetaLR_rankscore'][i]


### Round META LR score ###       
df_MAF['MetaLR_rankscore'] = pd.to_numeric(df_MAF['MetaLR_rankscore'], errors = 'coerce')
df_MAF['MetaLR_rankscore'] = df_MAF['MetaLR_rankscore'].round(decimals = 3)


### Select first META LR pred ###  
for i in range(len(df_MAF)):
    if isinstance(df_MAF['MetaLR_pred'][i], str):
        n = df_MAF['MetaLR_pred'][i].find(',')
        if n > 0:
            df_MAF['MetaLR_pred'][i] = df_MAF['MetaLR_pred'][i][0:n]
        elif n < 0:
            df_MAF['MetaLR_pred'][i] = df_MAF['MetaLR_pred'][i]


### Select first Reliability index Meta score ###
for i in range(len(df_MAF)):
    if isinstance(df_MAF['Reliability_index_Meta_score'][i], str):
        n = df_MAF['Reliability_index_Meta_score'][i].find(',')
        if n > 0:
            df_MAF['Reliability_index_Meta_score'][i] = df_MAF['Reliability_index_Meta_score'][i][0:n]
        elif n < 0:
            df_MAF['Reliability_index_Meta_score'][i] = df_MAF['Reliability_index_Meta_score'][i]


### Select first REVEL score ###
for i in range(len(df_MAF)):
    if isinstance(df_MAF['REVEL_score'][i], str):
        n = df_MAF['REVEL_score'][i].find(',')
        if n > 0:
            df_MAF['REVEL_score'][i] = df_MAF['REVEL_score'][i][0:n]
        elif n < 0:
            df_MAF['REVEL_score'][i] = df_MAF['REVEL_score'][i]


### Round REVEL score ### 
df_MAF['REVEL_score'] = pd.to_numeric(df_MAF['REVEL_score'], errors = 'coerce')
df_MAF['REVEL_score'] = df_MAF['REVEL_score'].round(decimals = 3)


### Select first REVEL rankscore ###
for i in range(len(df_MAF)):
    if isinstance(df_MAF['REVEL_rankscore'][i], str):
        n = df_MAF['REVEL_rankscore'][i].find(',')
        if n > 0:
            df_MAF['REVEL_rankscore'][i] = df_MAF['REVEL_rankscore'][i][0:n]
        elif n < 0:
            df_MAF['REVEL_rankscore'][i] = df_MAF['REVEL_rankscore'][i]


### Round REVEL rankscore ### 
df_MAF['REVEL_rankscore'] = pd.to_numeric(df_MAF['REVEL_rankscore'], errors = 'coerce')
df_MAF['REVEL_rankscore'] = df_MAF['REVEL_rankscore'].round(decimals = 3)


### Select first BayesDel AF score ###
for i in range(len(df_MAF)):
    if isinstance(df_MAF['BayesDel_addAF_score'][i], str):
        n = df_MAF['BayesDel_addAF_score'][i].find(',')
        if n > 0:
            df_MAF['BayesDel_addAF_score'][i] = df_MAF['BayesDel_addAF_score'][i][0:n]
        elif n < 0:
            df_MAF['BayesDel_addAF_score'][i] = df_MAF['BayesDel_addAF_score'][i]


### Round BayesDel AF score ### 
df_MAF['BayesDel_addAF_score'] = pd.to_numeric(df_MAF['BayesDel_addAF_score'], errors = 'coerce')
df_MAF['BayesDel_addAF_score'] = df_MAF['BayesDel_addAF_score'].round(decimals = 3)


### Select first BayesDel AF rankscore ### 
for i in range(len(df_MAF)):
    if isinstance(df_MAF['BayesDel_addAF_rankscore'][i], str):
        n = df_MAF['BayesDel_addAF_rankscore'][i].find(',')
        if n > 0:
            df_MAF['BayesDel_addAF_rankscore'][i] = df_MAF['BayesDel_addAF_rankscore'][i][0:n]
        elif n < 0:
            df_MAF['BayesDel_addAF_rankscore'][i] = df_MAF['BayesDel_addAF_rankscore'][i]


### Round BayesDel AF rankscore ### 
df_MAF['BayesDel_addAF_rankscore'] = pd.to_numeric(df_MAF['BayesDel_addAF_rankscore'], errors = 'coerce')
df_MAF['BayesDel_addAF_rankscore'] = df_MAF['BayesDel_addAF_rankscore'].round(decimals = 3)


### Select first BayesDel AF pred ### 
for i in range(len(df_MAF)):
    if isinstance(df_MAF['BayesDel_addAF_pred'][i], str):
        n = df_MAF['BayesDel_addAF_pred'][i].find(',')
        if n > 0:
            df_MAF['BayesDel_addAF_pred'][i] = df_MAF['BayesDel_addAF_pred'][i][0:n]
        elif n < 0:
            df_MAF['BayesDel_addAF_pred'][i] = df_MAF['BayesDel_addAF_pred'][i]


### Select first BayesDel no AF score ###
for i in range(len(df_MAF)):
    if isinstance(df_MAF['BayesDel_noAF_score'][i], str):
        n = df_MAF['BayesDel_noAF_score'][i].find(',')
        if n > 0:
            df_MAF['BayesDel_noAF_score'][i] = df_MAF['BayesDel_noAF_score'][i][0:n]
        elif n < 0:
            df_MAF['BayesDel_noAF_score'][i] = df_MAF['BayesDel_noAF_score'][i]


### Round BayesDel no AF score ###
df_MAF['BayesDel_noAF_score'] = pd.to_numeric(df_MAF['BayesDel_noAF_score'], errors = 'coerce')
df_MAF['BayesDel_noAF_score'] = df_MAF['BayesDel_noAF_score'].round(decimals = 3)


### Select first BayesDel no AF rankscore ### 
for i in range(len(df_MAF)):
    if isinstance(df_MAF['BayesDel_noAF_rankscore'][i], str):
        n = df_MAF['BayesDel_noAF_rankscore'][i].find(',')
        if n > 0:
            df_MAF['BayesDel_noAF_rankscore'][i] = df_MAF['BayesDel_noAF_rankscore'][i][0:n]
        elif n < 0:
            df_MAF['BayesDel_noAF_rankscore'][i] = df_MAF['BayesDel_noAF_rankscore'][i]


### Round BayesDel no AF rankscore ### 
df_MAF['BayesDel_noAF_rankscore'] = pd.to_numeric(df_MAF['BayesDel_noAF_rankscore'], errors = 'coerce')
df_MAF['BayesDel_noAF_rankscore'] = df_MAF['BayesDel_noAF_rankscore'].round(decimals = 3)


### Select first BayesDel no AF pred ### 
for i in range(len(df_MAF)):
    if isinstance(df_MAF['BayesDel_noAF_pred'][i], str):
        n = df_MAF['BayesDel_noAF_pred'][i].find(',')
        if n > 0:
            df_MAF['BayesDel_noAF_pred'][i] = df_MAF['BayesDel_noAF_pred'][i][0:n]
        elif n < 0:
            df_MAF['BayesDel_noAF_pred'][i] = df_MAF['BayesDel_noAF_pred'][i]


### Select first ClinPred score ### 
for i in range(len(df_MAF)):
    if isinstance(df_MAF['ClinPred_score'][i], str):
        n = df_MAF['ClinPred_score'][i].find(',')
        if n > 0:
            df_MAF['ClinPred_score'][i] = df_MAF['ClinPred_score'][i][0:n]
        elif n < 0:
            df_MAF['ClinPred_score'][i] = df_MAF['ClinPred_score'][i]


### Round ClinPred score ### 
df_MAF['ClinPred_score'] = pd.to_numeric(df_MAF['ClinPred_score'], errors = 'coerce')
df_MAF['ClinPred_score'] = df_MAF['ClinPred_score'].round(decimals = 3)


### Select first ClinPred rankscore ### 
for i in range(len(df_MAF)):
    if isinstance(df_MAF['ClinPred_rankscore'][i], str):
        n = df_MAF['ClinPred_rankscore'][i].find(',')
        if n > 0:
            df_MAF['ClinPred_rankscore'][i] = df_MAF['ClinPred_rankscore'][i][0:n]
        elif n < 0:
            df_MAF['ClinPred_rankscore'][i] = df_MAF['ClinPred_rankscore'][i]


### Round ClinPred rankscore ###  
df_MAF['ClinPred_rankscore'] = pd.to_numeric(df_MAF['ClinPred_rankscore'], errors = 'coerce')
df_MAF['ClinPred_rankscore'] = df_MAF['ClinPred_rankscore'].round(decimals = 3)


### Select first ClinPred pred ###  
for i in range(len(df_MAF)):
    if isinstance(df_MAF['ClinPred_pred'][i], str):
        n = df_MAF['ClinPred_pred'][i].find(',')
        if n > 0:
            df_MAF['ClinPred_pred'][i] = df_MAF['ClinPred_pred'][i][0:n]
        elif n < 0:
            df_MAF['ClinPred_pred'][i] = df_MAF['ClinPred_pred'][i]


### Select first DANN score ### 
for i in range(len(df_MAF)):
    if isinstance(df_MAF['DANN_score'][i], str):
        n = df_MAF['DANN_score'][i].find(',')
        if n > 0:
            df_MAF['DANN_score'][i] = df_MAF['DANN_score'][i][0:n]
        elif n < 0:
            df_MAF['DANN_score'][i] = df_MAF['DANN_score'][i]


### Round DANN score ### 
df_MAF['DANN_score'] = pd.to_numeric(df_MAF['DANN_score'], errors = 'coerce')
df_MAF['DANN_score'] = df_MAF['DANN_score'].round(decimals = 3)


### Select first DANN rankscore ###  
for i in range(len(df_MAF)):
    if isinstance(df_MAF['DANN_rankscore'][i], str):
        n = df_MAF['DANN_rankscore'][i].find(',')
        if n > 0:
            df_MAF['DANN_rankscore'][i] = df_MAF['DANN_rankscore'][i][0:n]
        elif n < 0:
            df_MAF['DANN_rankscore'][i] = df_MAF['DANN_rankscore'][i]


### Round DANN rankscore ### 
df_MAF['DANN_rankscore'] = pd.to_numeric(df_MAF['DANN_rankscore'], errors = 'coerce')
df_MAF['DANN_rankscore'] = df_MAF['DANN_rankscore'].round(decimals = 3)


### Convert Nan for variant type in string ### 
df_MAF['variant_type'] = df_MAF['variant_type'].replace(np.nan, '', regex = True)


### Select variants present in variant type file ### 
list_of_variants = []
i = 0
for i in range(len(df_MAF)):
    j = 0
    for j in range(len(variants_type)):
        if variants_type['variants'][j] in df_MAF['variant_type'][i]:
            list_of_variants.append(df_MAF['variant_type'][i])
        j += 1
    i += 1


### Get unique values from list of variants ###  
list_of_variants_unique = set(list_of_variants)


### Create new dataframe ###
df_variants = df_MAF.copy()
df_variants['variant_type_check'] = ""
i = 0


### Filter by type of variants ###
for i in range(len(df_variants["variant_type"])):
    if df_variants["variant_type"][i] in list_of_variants_unique:
        df_variants['variant_type_check'][i] = 'ok'
        if df_variants["variant_type"][i] == 'non_coding_transcript_exon_variant':
            df_variants['variant_type_check'][i] = ''


### Delete rows with no value in variant_type_check ### 
df_variants['variant_type_check'].replace('', np.nan, inplace = True)
df_variants.dropna(subset=['variant_type_check'], inplace = True)


### Delete columns in df_nm ### 
df_variants.drop(['variant_type_check'], inplace = True, axis = 1)


### Reset index ### 
df_variants = df_variants.reset_index(drop = True)


### Create columns with nm without version ### 
df_variants['temp_nm'] = ""
for i in range(len(df_variants)):
    n = df_variants['transcript_id'][i].find('.')
    df_variants['temp_nm'][i] = df_variants['transcript_id'][i][0:n]


### Merge df_variants with exon file ### 
df_variants_exons = df_variants.merge(exons, how = "left", left_on = "temp_nm", right_on = "RefSeq_mRNA_ID")


### Convert exon rank to integer ### 
df_variants_exons['exon_rank'] = df_variants_exons['exon_rank'].astype('Int64')
df_variants_exons['exon_total'] = df_variants_exons['exon_total'].astype('Int64')


### Convert exon rank and number in string ###  
df_variants_exons['exon_rank'] = df_variants_exons['exon_rank'].astype(str)
df_variants_exons['exon_total'] = df_variants_exons['exon_total'].astype(str)


### Create new column for exon rank ###  
df_variants_exons['exon'] = ""

for i in range(len(df_variants_exons)):
    z = df_variants_exons['exon_rank'][i] +"/" + df_variants_exons['exon_total'][i]
    df_variants_exons['exon'][i] = z


### Remove columns ### 
df_variants_exons = df_variants_exons.drop("gene_name", axis = 1)
df_variants_exons = df_variants_exons.drop("RefSeq_mRNA_ID", axis = 1)
df_variants_exons = df_variants_exons.drop("temp_nm", axis = 1)
df_variants_exons = df_variants_exons.drop("allele_depth", axis = 1)
df_variants_exons = df_variants_exons.drop("exon_rank", axis = 1)
df_variants_exons = df_variants_exons.drop("exon_total", axis = 1)


### Reorganize the order of columns ###  
df_variants_exons = df_variants_exons.reindex(['barcode', 'NIPP_TUMOSPEC', 'Centre', 'chromosome', 'position',
                                               'ref', 'alt', 'allele1_depth', 'allele2_depth',
                                               'sample_depth', 'allelic_frequency', 'variant_type',
                                               'gene_id', 'transcript_id', 'HGVS.c', 'HGVS.p',
                                               'exon', 'non_cancer_AF','non_cancer_AF_eas',
                                               'non_cancer_AF_asj', 'non_cancer_AF_afr',
                                               'non_cancer_AF_amr', 'non_cancer_AF_nfe',
                                               'non_cancer_AF_fin', 'non_cancer_AF_sas',
                                               'cadd_raw_v1.6', 'cadd_phred_v1.6', 'cadd_raw_v1.4',
                                               'cadd_phred_v1.4', 'VEST4_score', 'VEST4_rankscore', 
                                               'MetaSVM_score','MetaSVM_rankscore', 'MetaSVM_pred', 
                                               'MetaLR_score','MetaLR_rankscore', 'MetaLR_pred',
                                               'Reliability_index_Meta_score','REVEL_score', 'REVEL_rankscore', 
                                               'BayesDel_addAF_score','BayesDel_addAF_rankscore', 
                                               'BayesDel_addAF_pred','BayesDel_noAF_score', 'BayesDel_noAF_rankscore',
                                               'BayesDel_noAF_pred', 'ClinPred_score', 'ClinPred_rankscore',
                                               'ClinPred_pred', 'DANN_score', 'DANN_rankscore', 'clinvar_id',
                                               'clinvar_clnsig', 'clinvar_trait', 'clinvar_review',
                                               'clinvar_hgvs', 'clinvar_var_source', 'clinvar_MedGen_id',
                                               'clinvar_OMIM_id', 'clinvar_Orphanet_id'], axis = 1)


### Export final data frame in csv ### 
df_variants_exons.to_csv('df_variants.csv', index = False)


### Create new dataframe ### 
df_nm = df_variants_exons.copy()
df_nm['transcripts'] = ""
i = 0


### Formate NM column ###  
while i < len(df_nm)-1:
    j = 1
    x = str(df_nm['transcript_id'][i]) + ", " + str(df_nm['HGVS.c'][i]) + ", " + str(df_nm['HGVS.p'][i]) + ", " + str(df_nm['exon'][i])
    while i + j < len(df_nm) and df_nm['barcode'][i] == df_nm['barcode'][i+j] and df_nm['position'][i] == df_nm['position'][i+j] and df_nm['ref'][i] == df_nm['ref'][i+j] and df_nm['alt'][i] == df_nm['alt'][i+j]:
        x += ", " + str(df_nm['transcript_id'][i+j]) + ", " + str(df_nm['HGVS.c'][i+j]) + ", " + str(df_nm['HGVS.p'][i+j]) + ", " + str(df_nm['exon'][i+j])
        j += 1
    df_nm['transcripts'][i] = x
    i += j


### Delete rows with no value in new_nm ###  
df_nm['transcripts'].replace('', np.nan, inplace = True)
df_nm.dropna(subset=['transcripts'], inplace = True)


### Delete columns in df_nm ###  
df_nm.drop(['transcript_id', 'HGVS.c', 'HGVS.p', 'exon'], inplace = True, axis = 1)


### Reorganize columns in df_nm ###
df_nm = df_nm.reindex(['barcode', 'NIPP_TUMOSPEC', 'Centre', 'chromosome', 'position', 'ref', 'alt',
                       'allele1_depth', 'allele2_depth', 'sample_depth', 'allelic_frequency',
                       'variant_type', 'gene_id', 'transcripts', 'non_cancer_AF',
                       'non_cancer_AF_eas', 'non_cancer_AF_asj', 'non_cancer_AF_afr',
                       'non_cancer_AF_amr', 'non_cancer_AF_nfe', 'non_cancer_AF_fin',
                       'non_cancer_AF_sas', 'cadd_raw_v1.6', 'cadd_phred_v1.6',
                       'cadd_raw_v1.4', 'cadd_phred_v1.4', 'VEST4_score', 'VEST4_rankscore', 
                       'MetaSVM_score', 'MetaSVM_rankscore','MetaSVM_pred', 'MetaLR_score',
                       'MetaLR_rankscore', 'MetaLR_pred','Reliability_index_Meta_score', 
                       'REVEL_score', 'REVEL_rankscore','BayesDel_addAF_score',
                       'BayesDel_addAF_rankscore', 'BayesDel_addAF_pred', 'BayesDel_noAF_score', 
                       'BayesDel_noAF_rankscore', 'BayesDel_noAF_pred', 'ClinPred_score', 
                       'ClinPred_rankscore', 'ClinPred_pred', 'DANN_score', 'DANN_rankscore', 
                       'clinvar_id', 'clinvar_clnsig', 'clinvar_trait', 'clinvar_review', 
                       'clinvar_hgvs', 'clinvar_var_source', 'clinvar_MedGen_id',
                       'clinvar_OMIM_id', 'clinvar_Orphanet_id'], axis = 1)


### Reset index of df_nm ###
df_nm = df_nm.reset_index(drop = True)


### Export final data frame in csv ###
df_nm.to_csv('df_variants_formatted.csv', index = False)



### Number of variants for each step ###
# After concatenation
with open("number_of_variants.txt", "a") as file:
    file.write(f"Number of variants after concatenation : {df_concat.shape} \n")

# After filtering by allele frequency and read depth
with open("number_of_variants.txt", "a") as file:
    file.write(f"Number of variants after filtering by allele frequency and read depth : {df_AF_DP.shape} \n")

# After filtering by minor allelic frequency
with open("number_of_variants.txt", "a") as file:
    file.write(f"Number of variants after filtering by minor allelic frequency : {df_MAF.shape} \n")

# After filtering by types of variants
with open("number_of_variants.txt", "a") as file:
    file.write(f"Number of variants after filtering by types of variants : {df_variants.shape} \n")

# Number of unique variants
with open("number_of_variants.txt", "a") as file:
    file.write(f"Number of unique variants : {df_nm.shape} ")
