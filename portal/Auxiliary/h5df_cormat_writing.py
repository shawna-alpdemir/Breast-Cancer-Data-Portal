import numpy as np
import pandas as pd
import h5py
import time

# replace ~ with your directory
def Import_Correlation_Data():
    jo_protein_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/jo_pro_z.csv',index_col='Gene')
    jo_mrna_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/jo_mrna_z.csv',index_col='Gene')

    kr_protein_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/kr_pro_z.csv',index_col='Gene')
    kr_mrna_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/kr_rna_z.csv',index_col='Gene')

    me_protein_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/me_pro_z.csv',index_col='Gene')
    me_mrna_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/me_rna_z.csv',index_col='Gene')

    return jo_mrna_for_cor, kr_mrna_for_cor, me_mrna_for_cor, jo_protein_for_cor, kr_protein_for_cor, me_protein_for_cor

[jo_mrna_for_cor, kr_mrna_for_cor, me_mrna_for_cor, jo_protein_for_cor, kr_protein_for_cor, me_protein_for_cor] = Import_Correlation_Data()

perc = 30.0 # percentage of allowing NA
jo_min_count =  int(((100-perc)/100)*jo_protein_for_cor.shape[1] + 1)
jo_protein_for_cor = jo_protein_for_cor.dropna(axis=0,
                    thresh=jo_min_count)

kr_min_count =  int(((100-perc)/100)*kr_protein_for_cor.shape[1] + 1)
kr_protein_for_cor = kr_protein_for_cor.dropna(axis=0,
                    thresh=kr_min_count)

me_min_count =  int(((100-perc)/100)*jo_protein_for_cor.shape[1] + 1)
mod_df = jo_protein_for_cor.dropna(axis=0,
                    thresh=me_min_count)


start = time.time()

# create matrix
jo_protein_cormat = jo_protein_for_cor.transpose().corr()
kr_protein_cormat= kr_protein_for_cor.transpose().corr()
me_protein_cormat= me_protein_for_cor.transpose().corr()

# jo_mrna_cormat = jo_mrna_for_cor.transpose().corr()
# kr_mrna_cormat = kr_mrna_for_cor.transpose().corr()
# me_mrna_cormat = me_mrna_for_cor.transpose().corr()

# create hdf5 file
jo_protein_cormat_hdf5 = "/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/Clean HDF5 Cormat/jo_protein_cormat1.hdf5"
kr_protein_cormat_hdf5 = "/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/Clean HDF5 Cormat/kr_protein_cormat1.hdf5"
me_protein_cormat_hdf5 = "/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/Clean HDF5 Cormat/me_protein_cormat1.hdf5"
# jo_mrna_cormat_hdf5 = "/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/Clean HDF5 Cormat/jo_mrna_cormat.hdf5"
# kr_mrna_cormat_hdf5 = "/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/Clean HDF5 Cormat/kr_mrna_cormat.hdf5"
# me_mrna_cormat_hdf5 = "/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/Clean HDF5 Cormat/me_mrna_cormat.hdf5"


def write_hdf5(hdf5_FileName, DF_quant):
    """
    hdf5_FileName = the path to your hdf5 file
    DF_quant = the proteome/transcriptome dataset, with genes as index and columns as tumors
    DF_subtype = the group key dataset
    Subtype_col_num = integer, the column number where subtype names are stored
    """

    # change quantities matrix to np.array, tumors to list, gene index to list
    quantities = DF_quant.to_numpy()
    genes_in_column = DF_quant.columns.to_list()
    genes_in_index = DF_quant.index.to_list()

    # assemble dataframe
    DF = pd.DataFrame(data=quantities, index=genes_in_index, columns=genes_in_column)

    # write the DF into hdf5 file, containing each gene as a dataset and tumor list as dataset
    with h5py.File(hdf5_FileName, "w") as hdf5:
        for gene in genes_in_index:
            gene_quant = np.array(DF.loc[gene,:])
            hdf5.create_dataset(str(gene), data=gene_quant)

        hdf5.create_dataset('genes', data=genes_in_column)


# write the files
write_hdf5(jo_protein_cormat_hdf5, jo_protein_cormat)
write_hdf5(kr_protein_cormat_hdf5, kr_protein_cormat)
write_hdf5(me_protein_cormat_hdf5, me_protein_cormat)
# write_hdf5(jo_mrna_cormat_hdf5, jo_mrna_cormat)
# write_hdf5(kr_mrna_cormat_hdf5, kr_mrna_cormat)
# write_hdf5(me_mrna_cormat_hdf5, me_mrna_cormat)

end = time.time()
print(f"{end-start}")
