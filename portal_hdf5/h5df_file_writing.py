import h5py
import numpy as np
import pandas as pd

################################### Johansson
# create hdf5 files
JohanssonProteome = '/portal_hdf5/Data/JohanssonProteome.hdf5'
JohanssonTranscriptome = '/portal_hdf5/Data/JohanssonTranscriptome.hdf5'

# read in txt files
# jo_df_protein = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/ConvertedDateCleanData/jo_data_p.txt',
#                             sep='\t', index_col='Gene')
# jo_df_RNA = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/ConvertedDateCleanData/jo_data_m.txt',
#                             sep='\t', index_col='Gene')
# jo_df_subtype = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/ConvertedDateCleanData/jo_group_key.txt',
#                             sep='\t')

################################### Krug
# create hdf5 files
KrugProteome = '/portal_hdf5/Data/KrugProteome.hdf5'
KrugTranscriptome = '/portal_hdf5/Data/KrugTranscriptome.hdf5'

# read in txt files
# kr_df_protein = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/ConvertedDateCleanData/kr_data_p.txt',
#                             sep='\t', index_col='Gene')
# kr_df_RNA = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/ConvertedDateCleanData/kr_data_m.txt',
#                             sep='\t', index_col='Gene')
# kr_df_subtype = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/ConvertedDateCleanData/kr_group_key.txt',
#                             sep='\t')

################################### Mertins
# create hdf5 files
MertinsProteome = '/portal_hdf5/Data/MertinsProteome.hdf5'
MertinsTranscriptome = '/portal_hdf5/Data/MertinsTranscriptome.hdf5'

# read in txt files
# me_df_protein = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/ConvertedDateCleanData/me_data_p.txt',
#                             sep='\t', index_col='Gene')
# me_df_RNA = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/ConvertedDateCleanData/me_data_m.txt',
#                             sep='\t', index_col='Gene')
# me_df_subtype = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/ConvertedDateCleanData/me_group_key.txt',
#                             sep='\t')


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

    # aseemble dataframe
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
#
# write_hdf5(KrugProteome, kr_df_protein, kr_df_subtype, 2)
# write_hdf5(KrugTranscriptome, kr_df_RNA, kr_df_subtype, 2)
#
# write_hdf5(MertinsProteome, me_df_protein, me_df_subtype, 1)
# write_hdf5(MertinsTranscriptome, me_df_RNA, me_df_subtype, 1)

# test
with h5py.File(MertinsTranscriptome, "r") as hdf5:
    tumor_list = np.array(hdf5.get('tumors'))
    blank1_quants = np.array(hdf5.get('blank 1'))
    blank2_quants = np.array(hdf5.get('blank 2'))
    blank3_quants = np.array(hdf5.get('blank 3'))
    blank4_quants = np.array(hdf5.get('blank 4'))
    subtype_list = np.array(hdf5.get('subtypes'))
    gene_list = np.array(hdf5.get('genes'))
    # print(tumor_list)
    # print(len(subtype_list))
    print(type(blank1_quants))
    print(blank4_quants)

    print(gene_list)
    print(type(gene_list))
