import numpy as np
import pandas as pd
import h5py
import time

# replace ~ with your directory
def Import_Correlation_Data():
    jo_protein_for_cor = pd.read_csv('~/jo_data_p.txt',index_col='Gene', sep='\t')
    jo_mrna_for_cor = pd.read_csv('~/jo_data_m.txt',index_col='Gene', sep='\t')

    kr_protein_for_cor = pd.read_csv('~/kr_data_p.txt',index_col='Gene', sep='\t')
    kr_mrna_for_cor = pd.read_csv('~/kr_data_m.txt',index_col='Gene', sep='\t')

    me_protein_for_cor = pd.read_csv('~/me_data_p.txt',index_col='Gene', sep='\t')
    me_mrna_for_cor = pd.read_csv('~/me_data_m.txt',index_col='Gene', sep='\t')

    return jo_mrna_for_cor, kr_mrna_for_cor, me_mrna_for_cor, jo_protein_for_cor, kr_protein_for_cor, me_protein_for_cor

[jo_mrna_for_cor, kr_mrna_for_cor, me_mrna_for_cor, jo_protein_for_cor, kr_protein_for_cor, me_protein_for_cor] = Import_Correlation_Data()

start = time.time()

# create matrix
#jo_protein_cormat = jo_protein_for_cor.transpose().corr()
#kr_protein_cormat= kr_protein_for_cor.transpose().corr()
#me_protein_cormat= me_protein_for_cor.transpose().corr()

#jo_mrna_cormat = jo_mrna_for_cor.transpose().corr()
#kr_mrna_cormat = kr_mrna_for_cor.transpose().corr()
#me_mrna_cormat = me_mrna_for_cor.transpose().corr()

# create hdf5 file
#jo_protein_cormat_hdf5 = "~/jo_protein_cormat.hdf5"
#kr_protein_cormat_hdf5 = "~/kr_protein_cormat.hdf5"
#me_protein_cormat_hdf5 = "~/me_protein_cormat.hdf5"
#jo_mrna_cormat_hdf5 = "~/jo_mrna_cormat.hdf5"
#kr_mrna_cormat_hdf5 = "~/kr_mrna_cormat.hdf5"
#me_mrna_cormat_hdf5 = "~/me_mrna_cormat.hdf5"


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
#write_hdf5(jo_protein_cormat_hdf5, jo_protein_cormat)
#write_hdf5(kr_protein_cormat_hdf5, kr_protein_cormat)
#write_hdf5(me_protein_cormat_hdf5, me_protein_cormat)
#write_hdf5(jo_mrna_cormat_hdf5, jo_mrna_cormat)
#write_hdf5(kr_mrna_cormat_hdf5, kr_mrna_cormat)
#write_hdf5(me_mrna_cormat_hdf5, me_mrna_cormat)

end = time.time()
print(f"{end-start}")
