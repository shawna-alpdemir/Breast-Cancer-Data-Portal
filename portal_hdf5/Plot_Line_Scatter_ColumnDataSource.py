import h5py
import numpy as np
from bokeh.models import ColumnDataSource
from Import_Files import Import_HDF5
from pdb import set_trace

# function call
[JohanssonProteome, JohanssonTranscriptome, KrugProteome, KrugTranscriptome, MertinsProteome, MertinsTranscriptome] = Import_HDF5()

# dictionary key name
subtype_tuple = 'subtype, tumor'
subtype = 'subtype'
protein_data = 'protein_data'
mRNA_data = 'mRNA_data'
gene = 'gene'

########################################################################################################################
##################################################### Johansson ########################################################
def Johansson_CDS():
    """Create ColumnDataSource for subtype, subtype_tumor_tuple, mRNA, protein, gene list for plots with initial default genes"""

    with h5py.File(JohanssonProteome, "r") as a, h5py.File(JohanssonTranscriptome, "r") as b:
        tumor_list = np.array(a.get('tumors'))
        subtype_list = np.array(a.get('subtypes'))
        jo_pro_gene_list = np.array(a.get('genes'))
        jo_mrna_gene_list = np.array(b.get('genes'))

        tumor_list = [x.decode('utf-8') for x in tumor_list]  # decode bytes string
        subtype_list = [x.decode('utf-8') for x in subtype_list]
        jo_pro_gene_list = [x.decode('utf-8') for x in jo_pro_gene_list]
        jo_mrna_gene_list = [x.decode('utf-8') for x in jo_mrna_gene_list]

        jo_total_gene_list = jo_pro_gene_list + jo_mrna_gene_list
        jo_unique_gene_list = np.unique(jo_total_gene_list)
        jo_unique_gene_list = np.sort(jo_unique_gene_list) # return ndarray

        # create tuple containing subtype and tumor; for x_range in plotting
        johansson_subtype_tumor_tuple = tuple(zip(subtype_list, tumor_list))

        # Create the dictionaries for each line on the protein plot and fill them with values
        dict1 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS2')),
                 mRNA_data: np.array(b.get('NDUFS2')),
                 gene: np.repeat('NDUFS2', len(tumor_list))}
        dict2 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS3')),
                 mRNA_data: np.array(b.get('NDUFS3')),
                 gene: np.repeat('NDUFS3', len(tumor_list))}
        dict3 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS7')),
                 mRNA_data: np.array(b.get('NDUFS7')),
                 gene: np.repeat('NDUFS7', len(tumor_list))}
        dict4 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('blank 4')),
                 mRNA_data: np.array(b.get('blank 4')),
                 gene: np.repeat('blank 4', len(tumor_list))}
        dict5 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_protein_data': np.array(a.get('NDUFS2')),
                 'y_protein_data': np.array(a.get('NDUFS3'))}
        dict6 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_mRNA_data': np.array(b.get('NDUFS2')),
                 'y_mRNA_data': np.array(b.get('NDUFS3'))}
        johansson_cds = []
        for i in [dict1, dict2, dict3, dict4, dict5, dict6]:
            johansson_cds.append(ColumnDataSource(data=i))

    return johansson_cds, johansson_subtype_tumor_tuple, jo_unique_gene_list


########################################################################################################################
##################################################### Krug #############################################################
def Krug_CDS():
    """Create ColumnDataSource for mRNA and protein correlation line plot with initial default genes"""

    with h5py.File(KrugProteome, "r") as a, h5py.File(KrugTranscriptome, "r") as b:
        tumor_list = np.array(a.get('tumors'))
        subtype_list = np.array(a.get('subtypes'))
        kr_pro_gene_list = np.array(a.get('genes'))
        kr_mrna_gene_list = np.array(b.get('genes'))

        tumor_list = [x.decode('utf-8') for x in tumor_list]  # decode bytes string
        subtype_list = [x.decode('utf-8') for x in subtype_list]
        kr_pro_gene_list = [x.decode('utf-8') for x in kr_pro_gene_list]
        kr_mrna_gene_list = [x.decode('utf-8') for x in kr_mrna_gene_list]

        kr_total_gene_list = kr_pro_gene_list + kr_mrna_gene_list
        kr_unique_gene_list = np.unique(kr_total_gene_list)
        kr_unique_gene_list = np.sort(kr_unique_gene_list) # return ndarray

        krug_subtype_tumor_tuple = tuple(zip(subtype_list, tumor_list))
        # create tuple containing subtype and tumor; for x_range in plotting

        # Create the dictionaries for each line on the protein plot and fill them with values
        dict1 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS2')),
                 mRNA_data: np.array(b.get('NDUFS2')),
                 gene: np.repeat('NDUFS2', len(tumor_list))}
        dict2 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS3')),
                 mRNA_data: np.array(b.get('NDUFS3')),
                 gene: np.repeat('NDUFS3', len(tumor_list))}
        dict3 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS7')),
                 mRNA_data: np.array(b.get('NDUFS7')),
                 gene: np.repeat('NDUFS7', len(tumor_list))}
        dict4 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('blank 4')),
                 mRNA_data: np.array(b.get('blank 4')),
                 gene: np.repeat('blank 4', len(tumor_list))}
        dict5 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_protein_data': np.array(a.get('NDUFS2')),
                 'y_protein_data': np.array(a.get('NDUFS3'))}
        dict6 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_mRNA_data': np.array(b.get('NDUFS2')),
                 'y_mRNA_data': np.array(b.get('NDUFS3'))}
        krug_cds = []
        for i in [dict1, dict2, dict3, dict4, dict5, dict6]:
            krug_cds.append(ColumnDataSource(data=i))

    return krug_cds, krug_subtype_tumor_tuple, kr_unique_gene_list


########################################################################################################################
##################################################### Mertins ##########################################################
def Mertins_CDS():
    """Create ColumnDataSource for mRNA and protein correlation line plot with initial default genes"""

    with h5py.File(MertinsProteome, "r") as a, h5py.File(MertinsTranscriptome, "r") as b:
        tumor_list = np.array(a.get('tumors'))
        subtype_list = np.array(a.get('subtypes'))
        me_pro_gene_list = np.array(a.get('genes'))
        me_mrna_gene_list = np.array(b.get('genes'))

        tumor_list = [x.decode('utf-8') for x in tumor_list]  # decode bytes string
        subtype_list = [x.decode('utf-8') for x in subtype_list]
        me_pro_gene_list = [x.decode('utf-8') for x in me_pro_gene_list]
        me_mrna_gene_list = [x.decode('utf-8') for x in me_mrna_gene_list]

        me_total_gene_list = me_pro_gene_list + me_mrna_gene_list
        me_unique_gene_list = np.unique(me_total_gene_list)
        me_unique_gene_list = np.sort(me_unique_gene_list) # return ndarray

        mertins_subtype_tumor_tuple = tuple(zip(subtype_list, tumor_list))
        # create tuple containing subtype and tumor; for x_range in plotting

        # Create the dictionaries for each line on the protein plot and fill them with values
        dict1 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS2')),
                 mRNA_data: np.array(b.get('NDUFS2')),
                 gene: np.repeat('NDUFS2', len(tumor_list))}
        dict2 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS3')),
                 mRNA_data: np.array(b.get('NDUFS3')),
                 gene: np.repeat('NDUFS3', len(tumor_list))}
        dict3 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS7')),
                 mRNA_data: np.array(b.get('NDUFS7')),
                 gene: np.repeat('NDUFS7', len(tumor_list))}
        dict4 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('blank 4')),
                 mRNA_data: np.array(b.get('blank 4')),
                 gene: np.repeat('blank 4', len(tumor_list))}
        dict5 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_protein_data': np.array(a.get('NDUFS2')),
                 'y_protein_data': np.array(a.get('NDUFS3'))}
        dict6 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_mRNA_data': np.array(b.get('NDUFS2')),
                 'y_mRNA_data': np.array(b.get('NDUFS3'))}
        mertins_cds = []
        for i in [dict1, dict2, dict3, dict4, dict5, dict6]:
            mertins_cds.append(ColumnDataSource(data=i))

    return mertins_cds, mertins_subtype_tumor_tuple, me_unique_gene_list