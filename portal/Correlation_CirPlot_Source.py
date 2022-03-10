# import python packages
import numpy as np
from bokeh.models import ColumnDataSource

# import Import_data as im
# return jo_df_protein, jo_df_mrna, jo_genes, jo_subtypes, jo_subtype_colors

########################################################################################################################
##################################################### Johansson ########################################################

def jo_Cor2_Source(jo_df_protein, jo_df_mrna, INITIAL_GENE=None):
    """Dictionary source for mRNA-protein correlation plot with initial default genes"""
    # create an entry for a "blank" gene
    for i in ['blank 1','blank 2','blank 3','blank 4']:
        jo_df_protein[i] = np.nan
        jo_df_mrna[i] = np.nan

    # test this condition
    if INITIAL_GENE is None:
        INITIAL_GENE = ['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4']

    # Create the dictionaries for each line on the protein plot and fill them with the initial values
    jo_cor2_source = {'x1': jo_df_mrna.loc[:, INITIAL_GENE[0]], 'y1': jo_df_protein.loc[:, INITIAL_GENE[0]],
                      'type1': jo_df_protein['Pam50'],
                      'x2': jo_df_mrna.loc[:, INITIAL_GENE[1]], 'y2': jo_df_protein.loc[:, INITIAL_GENE[1]],
                      'type2': jo_df_protein['Pam50'],
                      'x3': jo_df_mrna.loc[:, INITIAL_GENE[2]], 'y3': jo_df_protein.loc[:, INITIAL_GENE[2]],
                      'type3': jo_df_protein['Pam50'],
                      'x4': jo_df_mrna.loc[:, INITIAL_GENE[3]], 'y4': jo_df_protein.loc[:, INITIAL_GENE[3]],
                      'type4': jo_df_protein['Pam50']}

    jo_cor2_source = ColumnDataSource(data=jo_cor2_source)

    return jo_cor2_source


########################################################################################################################
##################################################### Krug #############################################################

def kr_Cor2_Source(kr_df_protein, kr_df_mrna, INITIAL_GENE=None):
    """Dictionary source for mRNA-protein correlation plot with initial default genes"""
    # create an entry for a "blank" gene
    for i in ['blank 1','blank 2','blank 3','blank 4']:
        kr_df_protein[i] = np.nan
        kr_df_mrna[i] = np.nan

    # test this condition
    if INITIAL_GENE is None:
        INITIAL_GENE = ['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4']

    # Create the dictionaries for each line on the protein plot and fill them with the initial values
    kr_cor2_source = {'x1': kr_df_mrna.loc[:, INITIAL_GENE[0]], 'y1': kr_df_protein.loc[:, INITIAL_GENE[0]],
                      'type1': kr_df_protein['PAM50'],
                      'x2': kr_df_mrna.loc[:, INITIAL_GENE[1]], 'y2': kr_df_protein.loc[:, INITIAL_GENE[1]],
                      'type2': kr_df_protein['PAM50'],
                      'x3': kr_df_mrna.loc[:, INITIAL_GENE[2]], 'y3': kr_df_protein.loc[:, INITIAL_GENE[2]],
                      'type3': kr_df_protein['PAM50'],
                      'x4': kr_df_mrna.loc[:, INITIAL_GENE[3]], 'y4': kr_df_protein.loc[:, INITIAL_GENE[3]],
                      'type4': kr_df_protein['PAM50'], }

    kr_cor2_source = ColumnDataSource(data=kr_cor2_source)

    return kr_cor2_source


########################################################################################################################
##################################################### Mertins ##########################################################
def me_Cor2_Source(me_df_protein, me_df_mrna, INITIAL_GENE=None):
    """Dictionary source for mRNA-protein correlation plot with initial default genes"""
    # create an entry for a "blank" gene
    for i in ['blank 1','blank 2','blank 3','blank 4']:
        me_df_protein[i] = np.nan
        me_df_mrna[i] = np.nan

    # test this condition
    if INITIAL_GENE is None:
        INITIAL_GENE = ['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4']

    # Create the dictionaries for each line on the protein plot and fill them with the initial values
    me_cor2_source = {'x1': me_df_mrna.loc[:, INITIAL_GENE[0]], 'y1': me_df_protein.loc[:, INITIAL_GENE[0]],
                      'type1': me_df_protein['PAM50'],
                      'x2': me_df_mrna.loc[:, INITIAL_GENE[1]], 'y2': me_df_protein.loc[:, INITIAL_GENE[1]],
                      'type2': me_df_protein['PAM50'],
                      'x3': me_df_mrna.loc[:, INITIAL_GENE[2]], 'y3': me_df_protein.loc[:, INITIAL_GENE[2]],
                      'type3': me_df_protein['PAM50'],
                      'x4': me_df_mrna.loc[:, INITIAL_GENE[3]], 'y4': me_df_protein.loc[:, INITIAL_GENE[3]],
                      'type4': me_df_protein['PAM50'], }

    me_cor2_source = ColumnDataSource(data=me_cor2_source)

    return me_cor2_source

# # testing
# a=im.jo_ImportData()
# print(jo_init_Cor_Source(a[0], a[1]))
