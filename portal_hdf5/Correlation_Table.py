# this script is for creating the correlation table from protein and mRNA data
# with the hdf5 file open, try to get the data for the target gene,
# except target gene does not exist in the dataset, empty df will be returned instead.

import pandas as pd
import h5py
import numpy as np
import scipy.stats as sc
from Import_Files import Import_Cor_Matrix_HDF5
#import time

[jo_protein_cormat, jo_mrna_cormat, kr_protein_cormat, kr_mrna_cormat, me_protein_cormat, me_mrna_cormat] = Import_Cor_Matrix_HDF5()

# new method: use calculated correlation matrix HDF5 file to query gene, order and sort top 100, calculate p values respectively
def Get_Protein_Correlation_Table(gene_name):

    with h5py.File(jo_protein_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            DF.index.name = 'Gene'
            DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # only select observation 2 to 101, because first obs is itself.
            top_hundred = DF.iloc[1:101]

            # calculate p values
            TestStat = (top_hundred * (45 - 2) ** 0.5) / (1 - top_hundred ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat),df=45 - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF
            top_hundred.insert(1, 'p', Sig, True)
            top_hundred.reset_index(inplace=True)
            jo_Correlation_table_df = top_hundred

        except ValueError:
            jo_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)

    with h5py.File(kr_protein_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            DF.index.name = 'Gene'
            DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # only select observation 2 to 101, because first obs is itself.
            top_hundred = DF.iloc[1:101]

            # calculate p values
            TestStat = (top_hundred * (122 - 2) ** 0.5) / (1 - top_hundred ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat),df=122 - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF
            top_hundred.insert(1, 'p', Sig, True)
            top_hundred.reset_index(inplace=True)
            kr_Correlation_table_df = top_hundred

        except ValueError:
            kr_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)

    with h5py.File(me_protein_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            DF.index.name = 'Gene'
            DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # only select observation 2 to 101, because first obs is itself.
            top_hundred = DF.iloc[1:101]

            # calculate p values
            TestStat = (top_hundred * (77 - 2) ** 0.5) / (1 - top_hundred ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat),df=77 - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF
            top_hundred.insert(1, 'p', Sig, True)
            top_hundred.reset_index(inplace=True)
            me_Correlation_table_df = top_hundred

        except ValueError:
            me_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)

    return jo_Correlation_table_df, kr_Correlation_table_df, me_Correlation_table_df


def Get_mRNA_Correlation_Table(gene_name):
    with h5py.File(jo_mrna_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            DF.index.name = 'Gene'
            DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # only select observation 2 to 101, because first obs is itself.
            top_hundred = DF.iloc[1:101]

            # calculate p values
            TestStat = (top_hundred * (45 - 2) ** 0.5) / (1 - top_hundred ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat),
                          df=45 - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF
            top_hundred.insert(1, 'p', Sig, True)
            top_hundred.reset_index(inplace=True)
            jo_m_Correlation_table_df = top_hundred

        except ValueError:
            jo_m_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)

    with h5py.File(kr_mrna_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            DF.index.name = 'Gene'
            DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # only select observation 2 to 101, because first obs is itself.
            top_hundred = DF.iloc[1:101]

            # calculate p values
            TestStat = (top_hundred * (122 - 2) ** 0.5) / (1 - top_hundred ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat), df=122 - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF
            top_hundred.insert(1, 'p', Sig, True)
            top_hundred.reset_index(inplace=True)
            kr_m_Correlation_table_df = top_hundred

        except ValueError:
            kr_m_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)

    with h5py.File(me_mrna_cormat, "r") as a:

        try:
            # input gene
            gene = np.array(a.get(gene_name))

            # total genes
            total_gene_list = np.array(a.get('genes'))
            total_gene_list = [x.decode('utf-8') for x in total_gene_list]

            # assemble DF by putting correlation values in a column, with gene names in the index
            DF = pd.DataFrame(data=gene, index=total_gene_list, columns=['r'])
            DF.index.name = 'Gene'
            DF.sort_values(by='r', key=abs, inplace=True, ascending=False)

            # only select observation 2 to 101, because first obs is itself.
            top_hundred = DF.iloc[1:101]

            # calculate p values
            TestStat = (top_hundred * (77 - 2) ** 0.5) / (1 - top_hundred ** 2) ** 0.5  # calculate t score
            Sig = sc.t.sf(abs(TestStat), df=77 - 2) * 2  # using the survival function to compute one-sided p-value, then doubled it to make it two-side

            # assemble final DF
            top_hundred.insert(1, 'p', Sig, True)
            top_hundred.reset_index(inplace=True)
            me_m_Correlation_table_df = top_hundred

        except ValueError:
            me_m_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
                                                    "r": [0],
                                                    "p": [0]}, index=None)

    return jo_m_Correlation_table_df, kr_m_Correlation_table_df, me_m_Correlation_table_df