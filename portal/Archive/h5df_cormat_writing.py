# Archived Feb 4, 2023. Newest h5df writing file for correlation matrix is in jupyter notebook

import numpy as np
import pandas as pd
import h5py
import time

start = time.time()


# create hdf5 file
jo_pro_coef = "/PATH/jo_pro_coef.hdf5"
kr_pro_coef = "/PATH/kr_pro_coef.hdf5"
me_pro_coef = "/PATH/me_pro_coef.hdf5"

jo_rna_coef = "/PATH/jo_mrna_coef.hdf5"
kr_rna_coef = "/PATH/kr_mrna_coef.hdf5"
me_rna_coef = "/PATH/me_mrna_coef.hdf5"

jo_pro_pval = "/PATH/jo_pro_pval.hdf5"
kr_pro_pval = "/PATH/kr_pro_pval.hdf5"
me_pro_pval = "/PATH/me_pro_pval.hdf5"

jo_rna_pval = "/PATH/jo_mrna_pval.hdf5"
kr_rna_pval = "/PATH/kr_mrna_pval.hdf5"
me_rna_pval = "/PATH/me_mrna_pval.hdf5"

# read in the correlation matrices



def special_write_hdf5(hdf5_FileName, DF):
    """
    hdf5_FileName = the path to your hdf5 file
    DF_quant = the proteome/transcriptome dataset, with genes as index
    DF_subtype = the group key dataset
    Subtype_col_num = integer, the column number where subtype names are stored
    """

    # change quantities matrix to np.array, tumors to list, gene index to list
    # quantities = DF_quant.to_numpy()
    # genes_in_column = DF_quant.columns.to_list()
    # genes_in_index = DF_quant.index.to_list()

    # assemble dataframe
    # DF = pd.DataFrame(data=quantities, index=genes_in_index, columns=genes_in_column)

    start = time.time()

    genes_in_index = DF.index.to_list()
    genes_in_col = DF.columns.to_list()
    # write the DF into hdf5 file, containing each gene as a dataset and tumor list as dataset
    with h5py.File(hdf5_FileName, "w") as hdf5:
        for gene in genes_in_index:
            gene_quant = np.array(DF.loc[gene,:])
            hdf5.create_dataset(str(gene), data=gene_quant)

        hdf5.create_dataset('genes', data=genes_in_col)

    end = time.time()
    print(f"took {end - start} sec.")


# write the files
# write_hdf5(jo_pro_coef, jo_pro_coef_matrix)
# write_hdf5(kr_protein_cormat_hdf5, kr_protein_cormat)
# write_hdf5(me_protein_cormat_hdf5, me_protein_cormat)
# write_hdf5(jo_mrna_cormat_hdf5, jo_mrna_cormat)
# write_hdf5(kr_mrna_cormat_hdf5, kr_mrna_cormat)
# write_hdf5(me_mrna_cormat_hdf5, me_mrna_cormat)

end = time.time()
print(f"{end-start}")
