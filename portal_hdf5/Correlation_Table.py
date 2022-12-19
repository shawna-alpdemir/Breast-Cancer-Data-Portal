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

        except ValueError:
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

        except ValueError:
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

        except ValueError:
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
    protein_sum_df.drop([0],inplace=True) # drop the first gene because it is itself
    protein_sum_df.dropna(thresh=3, inplace=True) # save rows if NA is less than 2, thresh is 3 because you need gene with 2 p values (3 non NA values)

    protein_sum_df['p_jo'].mask(protein_sum_df['p_jo'] >= 0.01, 10000, inplace=True)
    protein_sum_df['p_jo'].mask(protein_sum_df['p_jo'] < 0.001, 1, inplace=True)
    protein_sum_df['p_jo'].mask(protein_sum_df['p_jo'] < 0.01, 10, inplace=True)
    protein_sum_df['p_jo'].mask(protein_sum_df['p_jo'] < 0.05, 100, inplace=True)
    protein_sum_df['p_jo'].mask(protein_sum_df['p_jo'] < 0.01, 1000, inplace=True)

    protein_sum_df['p_kr'].mask(protein_sum_df['p_kr'] >= 0.01, 10000, inplace=True)
    protein_sum_df['p_kr'].mask(protein_sum_df['p_kr'] < 0.001, 1, inplace=True)
    protein_sum_df['p_kr'].mask(protein_sum_df['p_kr'] < 0.01, 10, inplace=True)
    protein_sum_df['p_kr'].mask(protein_sum_df['p_kr'] < 0.05, 100, inplace=True)
    protein_sum_df['p_kr'].mask(protein_sum_df['p_kr'] < 0.01, 1000, inplace=True)

    protein_sum_df['p'].mask(protein_sum_df['p'] >= 0.01, 10000, inplace=True)
    protein_sum_df['p'].mask(protein_sum_df['p'] < 0.001, 1, inplace=True)
    protein_sum_df['p'].mask(protein_sum_df['p'] < 0.01, 10, inplace=True)
    protein_sum_df['p'].mask(protein_sum_df['p'] < 0.05, 100, inplace=True)
    protein_sum_df['p'].mask(protein_sum_df['p'] < 0.01, 1000, inplace=True)

    # sum the p values into a new column, and re-annotate with gene names, drop individual p values, rank by p-values
    protein_sum_df['p value'] = protein_sum_df.sum(numeric_only=True, axis=1)

    # category the p value levels
    for i in [2, 3, 12, 102, 1002, 10002]:
        protein_sum_df['p value'].mask(protein_sum_df['p value'] == i, '<0.001', inplace=True)
    for i in [11, 20, 21, 30, 111, 120, 1011, 1020, 10011, 10020]:
        protein_sum_df['p value'].mask(protein_sum_df['p value'] == i, '<0.01', inplace=True)
    for i in [101, 110, 200, 201, 210, 300, 1101, 1110, 1200, 10101, 10110, 10200]:
        protein_sum_df['p value'].mask(protein_sum_df['p value'] == i, '<0.05', inplace=True)
    for i in [1001, 1010, 1100, 2000, 2001, 2010, 2100, 3000, 11011, 11010, 11100, 12000]:
        protein_sum_df['p value'].mask(protein_sum_df['p value'] == i, '<0.1', inplace=True)
    for i in [10001, 10010, 10100, 11000, 20000, 20001, 20010, 20100, 21000, 30000]:
        protein_sum_df['p value'].mask(protein_sum_df['p value'] == i, '>0.1', inplace=True)

    protein_sum_df = protein_sum_df.join(annotate,on='Gene')
    protein_sum_df.set_index('Gene',inplace=True)
    protein_sum_df.drop(protein_sum_df.columns[[0,1,2]],axis=1,inplace=True)
    protein_sum_df['p value'] = pd.Categorical(protein_sum_df['p value'], ['<0.001', '<0.01', '<0.05', '<0.1', '>0.1'])
    protein_sum_df.sort_values(by='p value', inplace=True, ascending=True)

    protein_sum_df.to_csv("ERBB2_protein_summary_table.txt", sep='\t')
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

        except ValueError:
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

        except ValueError:
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

        except ValueError:
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
    mRNA_sum_df.drop([0],inplace=True) # drop the first gene because it is itself
    mRNA_sum_df.dropna(thresh=3, inplace=True) # save rows if NA is less than 2, thresh is 3 because you need gene with 2 p values (3 non NA values)

    mRNA_sum_df['p_jo'].mask(mRNA_sum_df['p_jo'] >= 0.01, 10000, inplace=True)
    mRNA_sum_df['p_jo'].mask(mRNA_sum_df['p_jo'] < 0.001, 1, inplace=True)
    mRNA_sum_df['p_jo'].mask(mRNA_sum_df['p_jo'] < 0.01, 10, inplace=True)
    mRNA_sum_df['p_jo'].mask(mRNA_sum_df['p_jo'] < 0.05, 100, inplace=True)
    mRNA_sum_df['p_jo'].mask(mRNA_sum_df['p_jo'] < 0.01, 1000, inplace=True)

    mRNA_sum_df['p_kr'].mask(mRNA_sum_df['p_kr'] >= 0.01, 10000, inplace=True)
    mRNA_sum_df['p_kr'].mask(mRNA_sum_df['p_kr'] < 0.001, 1, inplace=True)
    mRNA_sum_df['p_kr'].mask(mRNA_sum_df['p_kr'] < 0.01, 10, inplace=True)
    mRNA_sum_df['p_kr'].mask(mRNA_sum_df['p_kr'] < 0.05, 100, inplace=True)
    mRNA_sum_df['p_kr'].mask(mRNA_sum_df['p_kr'] < 0.01, 1000, inplace=True)

    mRNA_sum_df['p'].mask(mRNA_sum_df['p'] >= 0.01, 10000, inplace=True)
    mRNA_sum_df['p'].mask(mRNA_sum_df['p'] < 0.001, 1, inplace=True)
    mRNA_sum_df['p'].mask(mRNA_sum_df['p'] < 0.01, 10, inplace=True)
    mRNA_sum_df['p'].mask(mRNA_sum_df['p'] < 0.05, 100, inplace=True)
    mRNA_sum_df['p'].mask(mRNA_sum_df['p'] < 0.01, 1000, inplace=True)

    # sum the p values into a new column, and re-annotate with gene names, drop individual p values, rank by p-values
    mRNA_sum_df['p value'] = mRNA_sum_df.sum(numeric_only=True, axis=1)

    # category the p value levels
    for i in [2, 3, 12, 102, 1002, 10002]:
        mRNA_sum_df['p value'].mask(mRNA_sum_df['p value'] == i, '<0.001', inplace=True)
    for i in [11, 20, 21, 30, 111, 120, 1011, 1020, 10011, 10020]:
        mRNA_sum_df['p value'].mask(mRNA_sum_df['p value'] == i, '<0.01', inplace=True)
    for i in [101, 110, 200, 201, 210, 300, 1101, 1110, 1200, 10101, 10110, 10200]:
        mRNA_sum_df['p value'].mask(mRNA_sum_df['p value'] == i, '<0.05', inplace=True)
    for i in [1001, 1010, 1100, 2000, 2001, 2010, 2100, 3000, 11011, 11010, 11100, 12000]:
        mRNA_sum_df['p value'].mask(mRNA_sum_df['p value'] == i, '<0.1', inplace=True)
    for i in [10001, 10010, 10100, 11000, 20000, 20001, 20010, 20100, 21000, 30000]:
        mRNA_sum_df['p value'].mask(mRNA_sum_df['p value'] == i, '>0.1', inplace=True)

    # re-annotate with gene names, drop individual p values
    mRNA_sum_df = mRNA_sum_df.join(annotate,on='Gene')
    mRNA_sum_df.set_index('Gene',inplace=True)
    mRNA_sum_df.drop(mRNA_sum_df.columns[[0,1,2]],axis=1,inplace=True)
    mRNA_sum_df['p value'] = pd.Categorical(mRNA_sum_df['p value'], ['<0.001', '<0.01', '<0.05', '<0.1', '>0.1'])
    mRNA_sum_df.sort_values(by='p value', inplace=True, ascending=True)

    mRNA_sum_df.to_csv("ERBB2_mRNA_summary_table.txt", sep='\t')
    return jo_m_Correlation_table_df, kr_m_Correlation_table_df, me_m_Correlation_table_df, mRNA_sum_df

Get_Protein_Correlation_Table('ERBB2')
Get_mRNA_Correlation_Table('ERBB2')