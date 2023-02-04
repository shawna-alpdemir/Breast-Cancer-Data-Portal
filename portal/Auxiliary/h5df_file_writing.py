# create hdf5 files by writing each gene entry including blank genes into the file, and store gene name, tumor name,
# subtype information as a list. total 6 hdf5 files are created.

import h5py
import numpy as np
import pandas as pd
#import time

# replace ~ with your directory

################################### Johansson
# create hdf5 files
# JohanssonProteome = '/PATH/JohanssonProteome.hdf5'
# JohanssonTranscriptome = '/PATH/JohanssonTranscriptome.hdf5'

# read in txt files
# jo_df_protein = pd.read_csv('/PATH/jo_pro_z.csv', index_col='Gene')
# jo_df_RNA = pd.read_csv('/PATH/jo_mrna_z.csv', index_col='Gene')
# jo_df_subtype = pd.read_csv('/PATH/jo_group_key.txt', sep='\t')

################################### Krug
# create hdf5 files
# KrugProteome = '/PATH/KrugProteome.hdf5'
# KrugTranscriptome = '/PATH/KrugTranscriptome.hdf5'

# read in txt files
# kr_df_protein = pd.read_csv('/PATH/kr_pro_z.csv', index_col='Gene')
# kr_df_RNA = pd.read_csv('/PATH/kr_rna_z.csv',index_col='Gene')
# kr_df_subtype = pd.read_csv('/PATH/kr_group_key.txt', sep='\t')

################################### Mertins
# create hdf5 files
# MertinsProteome = '/PATH/MertinsProteome.hdf5'
# MertinsTranscriptome = '/PATH/MertinsTranscriptome.hdf5'

# read in txt files
# me_df_protein = pd.read_csv('/PATH/me_pro_z.csv', index_col='Gene')
# me_df_RNA = pd.read_csv('/PATH/me_rna_z.csv', index_col='Gene')
# me_df_subtype = pd.read_csv('/PATH/me_group_key.txt', sep='\t')


# function
def write_hdf5(hdf5_FileName, DF_quant, DF_subtype, Subtype_col_num):
    """
    hdf5_FileName = the path to your hdf5 file
    DF_quant = the proteome/transcriptome dataset, with genes as index and columns as tumors
    DF_subtype = the group key dataset
    Subtype_col_num = integer, the column number where subtype names are stored
    """

    # change quantities matrix to np.array, tumors to list, gene index to list
    quantities = DF_quant.to_numpy()
    tumors = DF_quant.columns.to_list()
    genes = DF_quant.index.to_list()
    subtypes = DF_subtype.iloc[:, Subtype_col_num].to_list()

    # blank entries
    # blank_entries = ['blank 1', 'blank 2', 'blank 3', 'blank 4']
    # genes.extend(blank_entries)
    #
    # blank_arrays = np.array([
    #     np.repeat(np.nan, len(tumors)),
    #     np.repeat(np.nan, len(tumors)),
    #     np.repeat(np.nan, len(tumors)),
    #     np.repeat(np.nan, len(tumors))])
    # quantities = np.vstack((quantities, blank_arrays))

    # assemble dataframe
    DF = pd.DataFrame(data=quantities, index=genes, columns=tumors)

    # write the DF into hdf5 file, containing each gene as a dataset and tumor list as dataset
    with h5py.File(hdf5_FileName, "w") as hdf5:
        for gene in genes:
            gene_quant = np.array(DF.loc[gene,:])
            hdf5.create_dataset(str(gene), data=gene_quant)

        hdf5.create_dataset('genes', data=genes)
        hdf5.create_dataset('tumors', data=tumors)
        hdf5.create_dataset('subtypes', data=subtypes)

# write the files
# write_hdf5(JohanssonProteome, jo_df_protein, jo_df_subtype, 2)
# write_hdf5(JohanssonTranscriptome, jo_df_RNA, jo_df_subtype, 2)

# write_hdf5(KrugProteome, kr_df_protein, kr_df_subtype, 2)
# write_hdf5(KrugTranscriptome, kr_df_RNA, kr_df_subtype, 2)

# write_hdf5(MertinsProteome, me_df_protein, me_df_subtype, 1)
# write_hdf5(MertinsTranscriptome, me_df_RNA, me_df_subtype, 1)
