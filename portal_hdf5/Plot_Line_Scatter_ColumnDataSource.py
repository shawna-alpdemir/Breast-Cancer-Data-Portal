# function that creates CDS for line plots and scatter plots

import h5py
import numpy as np
from bokeh.models import ColumnDataSource
from Import_Files import Import_HDF5
#from pdb import set_trace

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
        jo_dict1 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS2')),
                 mRNA_data: np.array(b.get('NDUFS2')),
                 gene: np.repeat('NDUFS2', len(tumor_list))}
        jo_dict2 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS3')),
                 mRNA_data: np.array(b.get('NDUFS3')),
                 gene: np.repeat('NDUFS3', len(tumor_list))}
        jo_dict3 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS7')),
                 mRNA_data: np.array(b.get('NDUFS7')),
                 gene: np.repeat('NDUFS7', len(tumor_list))}
        jo_dict4 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS8')),
                 mRNA_data: np.array(b.get('NDUFS8')),
                 gene: np.repeat('NDUFS8', len(tumor_list))}
        jo_dict5 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_protein_data': np.array(a.get('NDUFS2')),
                 'y_protein_data': np.array(a.get('NDUFS3'))}
        jo_dict6 = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_mRNA_data': np.array(b.get('NDUFS2')),
                 'y_mRNA_data': np.array(b.get('NDUFS3'))}
        johansson_cds = []
        for i in [jo_dict1, jo_dict2, jo_dict3, jo_dict4, jo_dict5, jo_dict6]:
            johansson_cds.append(ColumnDataSource(data=i))

    return johansson_cds, johansson_subtype_tumor_tuple, jo_unique_gene_list

def Johansson_Scatter_CDS():
    """Create ColumnDataSource for subtype, subtype_tumor_tuple, mRNA, protein, gene list for plots with initial default genes
    This one changes NA to zeros"""
    with h5py.File(JohanssonProteome, "r") as a, h5py.File(JohanssonTranscriptome, "r") as b:
        subtype_list = np.array(a.get('subtypes'))
        subtype_list = [x.decode('utf-8') for x in subtype_list]

        # Create the dictionaries for each line on the protein plot and fill them with values
        jo_dict1 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS2')),
                 mRNA_data: np.array(b.get('NDUFS2'))}
        jo_dict2 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS3')),
                 mRNA_data: np.array(b.get('NDUFS3'))}
        jo_dict3 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS7')),
                 mRNA_data: np.array(b.get('NDUFS7'))}
        jo_dict4 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS8')),
                 mRNA_data: np.array(b.get('NDUFS8'))}

        johansson_scatter_cds = []
        for i in [jo_dict1, jo_dict2, jo_dict3, jo_dict4]:
            johansson_scatter_cds.append(ColumnDataSource(data=i))

    return johansson_scatter_cds
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
        kr_dict1 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS2')),
                 mRNA_data: np.array(b.get('NDUFS2')),
                 gene: np.repeat('NDUFS2', len(tumor_list))}
        kr_dict2 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS3')),
                 mRNA_data: np.array(b.get('NDUFS3')),
                 gene: np.repeat('NDUFS3', len(tumor_list))}
        kr_dict3 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS7')),
                 mRNA_data: np.array(b.get('NDUFS7')),
                 gene: np.repeat('NDUFS7', len(tumor_list))}
        kr_dict4 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS8')),
                 mRNA_data: np.array(b.get('NDUFS8')),
                 gene: np.repeat('NDUFS8', len(tumor_list))}
        kr_dict5 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_protein_data': np.array(a.get('NDUFS2')),
                 'y_protein_data': np.array(a.get('NDUFS3'))}
        kr_dict6 = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_mRNA_data': np.array(b.get('NDUFS2')),
                 'y_mRNA_data': np.array(b.get('NDUFS3'))}
        krug_cds = []
        for i in [kr_dict1, kr_dict2, kr_dict3, kr_dict4, kr_dict5, kr_dict6]:
            krug_cds.append(ColumnDataSource(data=i))

    return krug_cds, krug_subtype_tumor_tuple, kr_unique_gene_list

def Krug_Scatter_CDS():
    """Create ColumnDataSource for mRNA and protein correlation line plot with initial default genes"""

    with h5py.File(KrugProteome, "r") as a, h5py.File(KrugTranscriptome, "r") as b:
        subtype_list = np.array(a.get('subtypes'))
        subtype_list = [x.decode('utf-8') for x in subtype_list]
        # create tuple containing subtype and tumor; for x_range in plotting

        # Create the dictionaries for each line on the protein plot and fill them with values
        kr_dict1 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS2')),
                 mRNA_data: np.array(b.get('NDUFS2'))}
        kr_dict2 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS3')),
                 mRNA_data: np.array(b.get('NDUFS3'))}
        kr_dict3 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS7')),
                 mRNA_data: np.array(b.get('NDUFS7'))}
        kr_dict4 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS8')),
                 mRNA_data: np.array(b.get('NDUFS8'))}

        krug_scatter_cds = []
        for i in [kr_dict1, kr_dict2, kr_dict3, kr_dict4]:
            krug_scatter_cds.append(ColumnDataSource(data=i))

    return krug_scatter_cds
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
        me_dict1 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS2')),
                 mRNA_data: np.array(b.get('NDUFS2')),
                 gene: np.repeat('NDUFS2', len(tumor_list))}
        me_dict2 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS3')),
                 mRNA_data: np.array(b.get('NDUFS3')),
                 gene: np.repeat('NDUFS3', len(tumor_list))}
        me_dict3 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS7')),
                 mRNA_data: np.array(b.get('NDUFS7')),
                 gene: np.repeat('NDUFS7', len(tumor_list))}
        me_dict4 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS8')),
                 mRNA_data: np.array(b.get('NDUFS8')),
                 gene: np.repeat('NDUFS8', len(tumor_list))}
        me_dict5 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_protein_data': np.array(a.get('NDUFS2')),
                 'y_protein_data': np.array(a.get('NDUFS3'))}
        me_dict6 = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: subtype_list,
                 'x_mRNA_data': np.array(b.get('NDUFS2')),
                 'y_mRNA_data': np.array(b.get('NDUFS3'))}
        mertins_cds = []
        for i in [me_dict1, me_dict2, me_dict3, me_dict4, me_dict5, me_dict6]:
            mertins_cds.append(ColumnDataSource(data=i))

    return mertins_cds, mertins_subtype_tumor_tuple, me_unique_gene_list

def Mertins_Scatter_CDS():
    """Create ColumnDataSource for mRNA and protein correlation line plot with initial default genes"""

    with h5py.File(MertinsProteome, "r") as a, h5py.File(MertinsTranscriptome, "r") as b:
        subtype_list = np.array(a.get('subtypes'))
        subtype_list = [x.decode('utf-8') for x in subtype_list]

        # Create the dictionaries for each line on the protein plot and fill them with values
        me_dict1 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS2')),
                 mRNA_data: np.array(b.get('NDUFS2'))}
        me_dict2 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS3')),
                 mRNA_data: np.array(b.get('NDUFS3'))}
        me_dict3 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS7')),
                 mRNA_data: np.array(b.get('NDUFS7'))}
        me_dict4 = {subtype: subtype_list,
                 protein_data: np.array(a.get('NDUFS8')),
                 mRNA_data: np.array(b.get('NDUFS8'))}

        mertins_scatter_cds = []
        for i in [me_dict1, me_dict2, me_dict3, me_dict4]:
            mertins_scatter_cds.append(ColumnDataSource(data=i))

    return mertins_scatter_cds