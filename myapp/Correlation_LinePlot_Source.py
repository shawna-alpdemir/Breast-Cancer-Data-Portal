import numpy as np
from bokeh.models import ColumnDataSource
# import Import_data as im

subtype = 'subtype'
protein_data = 'protein_data'
mRNA_data = 'mRNA_data'
gene = 'gene'

########################################################################################################################
##################################################### Johansson ########################################################

def jo_Cor1_Source(jo_df_protein, jo_df_mrna, INITIAL_GENE=None):
    """Dictionary source for mRNA and protein correlation line plot with initial default genes"""
    # create an entry for a "blank" gene
    for i in ['blank 1','blank 2','blank 3','blank 4']:
        jo_df_protein[i] = np.nan
        jo_df_mrna[i] = np.nan

    # test these conditions
    if INITIAL_GENE is None:
        INITIAL_GENE = ['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4']

    # Create the dictionaries for each line on the protein plot and fill them with values
    jo_cor1_source_dict1 = {subtype: tuple(zip(jo_df_protein['Pam50'], jo_df_protein.index)),
                            protein_data: jo_df_protein.loc[:, INITIAL_GENE[0]], mRNA_data: jo_df_mrna.loc[:, INITIAL_GENE[0]],
                            gene: np.repeat('NDUFS2', len(jo_df_protein['Pam50']))}
    jo_cor1_source_dict2 = {subtype: tuple(zip(jo_df_protein['Pam50'], jo_df_protein.index)),
                            protein_data: jo_df_protein.loc[:, INITIAL_GENE[1]], mRNA_data: jo_df_mrna.loc[:, INITIAL_GENE[1]],
                            gene: np.repeat('NDUFS3', len(jo_df_protein['Pam50']))}
    jo_cor1_source_dict3 = {subtype: tuple(zip(jo_df_protein['Pam50'], jo_df_protein.index)),
                            protein_data: jo_df_protein.loc[:, INITIAL_GENE[2]], mRNA_data: jo_df_mrna.loc[:, INITIAL_GENE[2]],
                            gene: np.repeat('NDUFS7', len(jo_df_protein['Pam50']))}
    jo_cor1_source_dict4 = {subtype: tuple(zip(jo_df_protein['Pam50'], jo_df_protein.index)),
                            protein_data: jo_df_protein.loc[:, INITIAL_GENE[3]], mRNA_data: jo_df_mrna.loc[:, INITIAL_GENE[3]],
                            gene: np.repeat('blank 4', len(jo_df_protein['Pam50']))}

    jo_SourceList = []
    for i in [jo_cor1_source_dict1, jo_cor1_source_dict2, jo_cor1_source_dict3, jo_cor1_source_dict4]:
        jo_SourceList.append(ColumnDataSource(data=i))

    return jo_SourceList


########################################################################################################################
##################################################### Krug #############################################################

def kr_Cor1_Source(kr_df_protein, kr_df_mrna, INITIAL_GENE=None):
    """Dictionary source for mRNA and protein correlation line plot with initial default genes"""
    # create an entry for a "blank" gene
    for i in ['blank 1','blank 2','blank 3','blank 4']:
        kr_df_protein[i] = np.nan
        kr_df_mrna[i] = np.nan

    # test this condition
    if INITIAL_GENE is None:
        INITIAL_GENE = ['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4']

    # Create the dictionaries for each line on the protein plot and fill them with values
    kr_cor1_source_dict1 = {subtype: tuple(zip(kr_df_protein['PAM50'], kr_df_protein.index)),
                            protein_data: kr_df_protein.loc[:, INITIAL_GENE[0]], mRNA_data: kr_df_mrna.loc[:, INITIAL_GENE[0]],
                            gene: np.repeat('NDUFS2', len(kr_df_protein['PAM50']))}
    kr_cor1_source_dict2 = {subtype: tuple(zip(kr_df_protein['PAM50'], kr_df_protein.index)),
                            protein_data: kr_df_protein.loc[:, INITIAL_GENE[1]], mRNA_data: kr_df_mrna.loc[:, INITIAL_GENE[1]],
                            gene: np.repeat('NDUFS3', len(kr_df_protein['PAM50']))}
    kr_cor1_source_dict3 = {subtype: tuple(zip(kr_df_protein['PAM50'], kr_df_protein.index)),
                            protein_data: kr_df_protein.loc[:, INITIAL_GENE[2]], mRNA_data: kr_df_mrna.loc[:, INITIAL_GENE[2]],
                            gene: np.repeat('NDUFS7', len(kr_df_protein['PAM50']))}
    kr_cor1_source_dict4 = {subtype: tuple(zip(kr_df_protein['PAM50'], kr_df_protein.index)),
                            protein_data: kr_df_protein.loc[:, INITIAL_GENE[3]], mRNA_data: kr_df_mrna.loc[:, INITIAL_GENE[3]],
                            gene: np.repeat('blank 4', len(kr_df_protein['PAM50']))}

    kr_SourceList = []
    for i in [kr_cor1_source_dict1, kr_cor1_source_dict2, kr_cor1_source_dict3, kr_cor1_source_dict4]:
        kr_SourceList.append(ColumnDataSource(data=i))

    return kr_SourceList


########################################################################################################################
##################################################### Mertins ##########################################################

def me_Cor1_Source(me_df_protein, me_df_mrna, INITIAL_GENE=None):
    """Dictionary source for mRNA and protein correlation line plot with initial default genes"""
    # create an entry for a "blank" gene
    for i in ['blank 1','blank 2','blank 3','blank 4']:
        me_df_protein[i] = np.nan
        me_df_mrna[i] = np.nan

    # test this condition
    if INITIAL_GENE is None:
        INITIAL_GENE = ['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4']

    # Create the dictionaries for each line on the protein plot and fill them with values
    me_cor1_source_dict1 = {subtype: tuple(zip(me_df_protein['PAM50'], me_df_protein.index)),
                            protein_data: me_df_protein.loc[:, INITIAL_GENE[0]], mRNA_data: me_df_mrna.loc[:, INITIAL_GENE[0]],
                            gene: np.repeat('NDUFS2', len(me_df_protein['PAM50']))}
    me_cor1_source_dict2 = {subtype: tuple(zip(me_df_protein['PAM50'], me_df_protein.index)),
                            protein_data: me_df_protein.loc[:, INITIAL_GENE[1]], mRNA_data: me_df_mrna.loc[:, INITIAL_GENE[1]],
                            gene: np.repeat('NDUFS3', len(me_df_protein['PAM50']))}
    me_cor1_source_dict3 = {subtype: tuple(zip(me_df_protein['PAM50'], me_df_protein.index)),
                            protein_data: me_df_protein.loc[:, INITIAL_GENE[2]], mRNA_data: me_df_mrna.loc[:, INITIAL_GENE[2]],
                            gene: np.repeat('NDUFS7', len(me_df_protein['PAM50']))}
    me_cor1_source_dict4 = {subtype: tuple(zip(me_df_protein['PAM50'], me_df_protein.index)),
                            protein_data: me_df_protein.loc[:, INITIAL_GENE[3]], mRNA_data: me_df_mrna.loc[:, INITIAL_GENE[3]],
                            gene: np.repeat('blank 4', len(me_df_protein['PAM50']))}

    me_SourceList = []
    for i in [me_cor1_source_dict1, me_cor1_source_dict2, me_cor1_source_dict3, me_cor1_source_dict4]:
        me_SourceList.append(ColumnDataSource(data=i))

    return me_SourceList

# # testing
# a = im.jo_ImportData()
# print(jo_Cor1_Source(a[0],a[1])[3])
