# function that create total gene list, all protein genes and mRNA genes are combined

# import
import numpy as np
from Plot_Line_Scatter_ColumnDataSource import Johansson_CDS, Krug_CDS, Mertins_CDS

# function calls
[johansson_cds, johansson_subtype_tumor_tuple, jo_unique_gene_list] = Johansson_CDS()
[krug_cds, krug_subtype_tumor_tuple, kr_unique_gene_list] = Krug_CDS()
[mertins_cds, mertins_subtype_tumor_tuple, me_unique_gene_list] = Mertins_CDS()

def Gene_List():
    all_unique_genes = np.concatenate([jo_unique_gene_list,kr_unique_gene_list,me_unique_gene_list])
    all_unique_genes = np.unique(all_unique_genes)
    all_unique_genes = np.sort(all_unique_genes)
    all_unique_genes = all_unique_genes.tolist()
    return all_unique_genes # 37997 entries, list type
