# this script is for creating the correlation table from protein and mRNA data
# with the hdf5 file open, try to get the data for the target gene,
# except target gene does not exist in the dataset, empty df will be returned instead.

import pandas as pd
import h5py
import numpy as np
import scipy.stats as sc
import pdb
from Import_Files import Import_Cor_Matrix_HDF5

jo_patients_num = 45
kr_patients_num = 122
me_patients_num = 77
[jo_protein_cormat, jo_mrna_cormat, kr_protein_cormat, kr_mrna_cormat, me_protein_cormat, me_mrna_cormat, annotate] = Import_Cor_Matrix_HDF5()

# new method: use calculated correlation matrix HDF5 file to query gene, order and sort top 1000, calculate p values respectively
def Get_Protein_Correlation_Table(gene_name):

    with h5py.File(jo_protein_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            jo_pro_DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            jo_pro_DF.index.name = 'Gene'
            jo_pro_DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # calculate p values
            TestStat = (jo_pro_DF * (jo_patients_num - 2) ** 0.5) / (1 - jo_pro_DF ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat),df=jo_patients_num - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF, X_pro_DF is for summary table assembly
            jo_pro_DF.insert(1, 'p', Sig, True)
            jo_pro_DF.reset_index(inplace=True)
            jo_Correlation_table_df = jo_pro_DF.iloc[1:1001]
            jo_Correlation_table_df = jo_Correlation_table_df.join(annotate, on='Gene')

        except ValueError or KeyError:
            jo_pro_DF = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)
            jo_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0],
                                                    "Name": [0]}, index=None)

    with h5py.File(kr_protein_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            kr_pro_DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            kr_pro_DF.index.name = 'Gene'
            kr_pro_DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # calculate p values
            TestStat = (kr_pro_DF * (kr_patients_num - 2) ** 0.5) / (1 - kr_pro_DF ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat),df=kr_patients_num - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF, X_pro_DF is for summary table assembly
            kr_pro_DF.insert(1, 'p', Sig, True)
            kr_pro_DF.reset_index(inplace=True)
            kr_Correlation_table_df = kr_pro_DF.iloc[1:1001]
            kr_Correlation_table_df = kr_Correlation_table_df.join(annotate, on='Gene')

        except ValueError or KeyError:
            kr_pro_DF = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)
            kr_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0],
                                                    "Name": [0]}, index=None)

    with h5py.File(me_protein_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            me_pro_DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            me_pro_DF.index.name = 'Gene'
            me_pro_DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # calculate p values
            TestStat = (me_pro_DF * (me_patients_num - 2) ** 0.5) / (1 - me_pro_DF ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat),df=me_patients_num - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF, X_pro_DF is for summary table assembly
            me_pro_DF.insert(1, 'p', Sig, True)
            me_pro_DF.reset_index(inplace=True)
            me_Correlation_table_df = me_pro_DF.iloc[1:1001]
            me_Correlation_table_df = me_Correlation_table_df.join(annotate, on='Gene')

        except ValueError or KeyError:
            me_pro_DF = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)
            me_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0],
                                                    "Name":[0]}, index=None)

    # remove coefficient column (col 1)
    jo_pro_DF.drop(jo_pro_DF.columns[[1]],axis=1,inplace=True)
    kr_pro_DF.drop(kr_pro_DF.columns[[1]],axis=1,inplace=True)
    me_pro_DF.drop(me_pro_DF.columns[[1]],axis=1,inplace=True)

    # merge dataframes, drop genes with only 1 p value
    jokr_pro_df = jo_pro_DF.merge(kr_pro_DF, how='outer', on='Gene', suffixes=['_jo','_kr'])
    protein_sum_df = jokr_pro_df.merge(me_pro_DF, how='outer', on='Gene')

    try:
        protein_sum_df.set_index('Gene', inplace=True)
        protein_sum_df.drop(labels=gene_name,inplace=True) # drop the first gene because it is itself
        protein_sum_df.reset_index(inplace=True)
    except KeyError: # skip the dropping if the gene has no data
        pass

    protein_sum_df.dropna(thresh=3, inplace=True) # save rows if NA is less than 2, thresh is 3 because you need gene with 2 p values (3 non NA values)

    # subset gene with 3 p values into df_3study, have 2 p value into df_2study
    protein_sum_df_3study = protein_sum_df.dropna(how='any')
    protein_sum_df_2study = pd.DataFrame(protein_sum_df.loc[protein_sum_df.isnull().any(axis=1)])

    # replace max p value of each gene with NaN
    protein_sum_df_3study = protein_sum_df_3study.mask(
        protein_sum_df_3study.eq(protein_sum_df_3study.max(axis=1, numeric_only=True), axis=0),
        np.nan,
        axis=0)

    # average the p values into a new column, merge the df back together
    protein_sum_df_3study['p avg'] = protein_sum_df_3study.mean(numeric_only=True, axis=1)
    protein_sum_df_2study['p avg'] = protein_sum_df_2study.mean(numeric_only=True, axis=1)
    protein_sum_df = protein_sum_df_2study.merge(protein_sum_df_3study, how='outer').sort_values('p avg')

    # re-annotate with gene names, drop individual p values, rank by p-values
    protein_sum_df = protein_sum_df.join(annotate,on='Gene')
    protein_sum_df.set_index('Gene',inplace=True)
    protein_sum_df.drop(protein_sum_df.columns[[0,1,2]],axis=1,inplace=True)

    return jo_Correlation_table_df, kr_Correlation_table_df, me_Correlation_table_df, protein_sum_df

def Get_mRNA_Correlation_Table(gene_name):

    with h5py.File(jo_mrna_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            jo_mRNA_DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            jo_mRNA_DF.index.name = 'Gene'
            jo_mRNA_DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # calculate p values
            TestStat = (jo_mRNA_DF * (jo_patients_num - 2) ** 0.5) / (1 - jo_mRNA_DF ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat), df=jo_patients_num - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF
            jo_mRNA_DF.insert(1, 'p', Sig, True)
            jo_mRNA_DF.reset_index(inplace=True)
            jo_m_Correlation_table_df = jo_mRNA_DF.iloc[1:1001]
            jo_m_Correlation_table_df = jo_m_Correlation_table_df.join(annotate, on='Gene')

        except ValueError or KeyError:
            jo_mRNA_DF = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)
            jo_m_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0],
                                                    "Name": [0]}, index=None)

    with h5py.File(kr_mrna_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            kr_mRNA_DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            kr_mRNA_DF.index.name = 'Gene'
            kr_mRNA_DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # calculate p values
            TestStat = (kr_mRNA_DF * (kr_patients_num - 2) ** 0.5) / (1 - kr_mRNA_DF ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat), df=kr_patients_num - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF
            kr_mRNA_DF.insert(1, 'p', Sig, True)
            kr_mRNA_DF.reset_index(inplace=True)
            kr_m_Correlation_table_df = kr_mRNA_DF.iloc[1:1001]
            kr_m_Correlation_table_df = kr_m_Correlation_table_df.join(annotate, on='Gene')

        except ValueError or KeyError:
            kr_mRNA_DF = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)
            kr_m_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0],
                                                    "Name": [0]}, index=None)

    with h5py.File(me_mrna_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            me_mRNA_DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            me_mRNA_DF.index.name = 'Gene'
            me_mRNA_DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # calculate p values
            TestStat = (me_mRNA_DF * (me_patients_num - 2) ** 0.5) / (1 - me_mRNA_DF ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat), df=me_patients_num - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF
            me_mRNA_DF.insert(1, 'p', Sig, True)
            me_mRNA_DF.reset_index(inplace=True)
            me_m_Correlation_table_df = me_mRNA_DF.iloc[1:1001]
            me_m_Correlation_table_df = me_m_Correlation_table_df.join(annotate, on='Gene')

        except ValueError or KeyError:
            me_mRNA_DF = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)
            me_m_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0],
                                                    "Name": [0]}, index=None)

    # remove coefficient columns
    jo_mRNA_DF.drop(jo_mRNA_DF.columns[[1]],axis=1,inplace=True)
    kr_mRNA_DF.drop(kr_mRNA_DF.columns[[1]],axis=1,inplace=True)
    me_mRNA_DF.drop(me_mRNA_DF.columns[[1]],axis=1,inplace=True)

    # merge dataframes, drop genes with only 1 p value
    jokr_mRNA_df = jo_mRNA_DF.merge(kr_mRNA_DF, how='outer', on='Gene', suffixes=['_jo','_kr'])
    mRNA_sum_df = jokr_mRNA_df.merge(me_mRNA_DF, how='outer', on='Gene')

    try:
        mRNA_sum_df.set_index('Gene', inplace=True)
        mRNA_sum_df.drop(labels=gene_name,inplace=True) # drop the first gene because it is itself
        mRNA_sum_df.reset_index(inplace=True)
    except KeyError: # skip the dropping if the gene has no data
        pass

    mRNA_sum_df.dropna(thresh=3, inplace=True) # save rows if NA is less than 2, thresh is 3 because you need gene with 2 p values (3 non NA values)

    # subset gene with 3 p values into df_3study, have 2 p value into df_2study
    mRNA_sum_df_3study = mRNA_sum_df.dropna(how='any')
    mRNA_sum_df_2study = pd.DataFrame(mRNA_sum_df.loc[mRNA_sum_df.isnull().any(axis=1)])

    # replace max p value of each gene with NaN
    mRNA_sum_df_3study = mRNA_sum_df_3study.mask(
        mRNA_sum_df_3study.eq(mRNA_sum_df_3study.max(axis=1, numeric_only=True), axis=0),
        np.nan,
        axis=0)

    # average the p values into a new column, merge the df back together
    mRNA_sum_df_3study['p avg'] = mRNA_sum_df_3study.mean(numeric_only=True, axis=1)
    mRNA_sum_df_2study['p avg'] = mRNA_sum_df_2study.mean(numeric_only=True, axis=1)
    mRNA_sum_df = mRNA_sum_df_2study.merge(mRNA_sum_df_3study, how='outer').sort_values('p avg')

    # re-annotate with gene names, drop individual p values, rank by p-values
    mRNA_sum_df = mRNA_sum_df.join(annotate,on='Gene')
    mRNA_sum_df.set_index('Gene',inplace=True)
    mRNA_sum_df.drop(mRNA_sum_df.columns[[0,1,2]],axis=1,inplace=True)

    return jo_m_Correlation_table_df, kr_m_Correlation_table_df, me_m_Correlation_table_df, mRNA_sum_df