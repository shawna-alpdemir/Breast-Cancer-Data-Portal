import numpy as np
import pandas as pd

########################################################################################################################
########################################## Import Data from Johansson ##################################################
def jo_ImportData():
    """
    Import protein, mRNA, metadata from Johansson et al
    """
    # annotation file
    jo_annotation = pd.read_csv('myapp/Data/Johansson/jo_group_key.txt', sep='\t', header='infer')
    jo_annotation = jo_annotation.set_index('Patient')  # set annotation file index as sample name

    # protein file
    jo_df_protein = pd.read_csv('myapp/Data/Johansson/jo_data_p.txt', sep='\t')

    # mRNA file
    jo_df_mrna = pd.read_csv('myapp/Data/Johansson/jo_data_m.txt', sep='\t')

    # unique gene list
    jo_genes = pd.unique(
        list(jo_df_protein.Gene) + list(jo_df_mrna.Gene))  # combine the genes, this will produce an array
    jo_genes = np.ndarray.tolist(jo_genes)  # convert the array to a list of genes that are not duplicated.

    # transpose dataset
    jo_df_protein = jo_df_protein.rename(columns={'Gene': 'Patient'}).set_index(
        'Patient').T  # rename the column for transposed dataset
    jo_df_mrna = jo_df_mrna.rename(columns={'Gene': 'Patient'}).set_index(
        'Patient').T  # rename the column for transposed dataset

    # subtype information
    jo_subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']  # list of subtypes
    jo_subtype_colors = ['#E31A1C', '#FB9A99', '#1F78B4', '#A6CEE3',
                         '#33A02C']  # list of corresponding colors for subtypes

    jo_df_protein = jo_df_protein.join(jo_annotation.iloc[:, [1]])  # merge subtypes info onto protein dataset
    jo_df_mrna = jo_df_mrna.join(jo_annotation.iloc[:, [1]])  # merge subtypes info onto mRNA dataset

    # sort data by subtype
    jo_df_protein = jo_df_protein.sort_index().sort_values(by='Pam50', kind='mergesort', ascending=True)
    jo_df_mrna = jo_df_mrna.sort_index().sort_values(by='Pam50', kind='mergesort', ascending=True)

    return jo_df_protein, jo_df_mrna, jo_genes, jo_subtypes, jo_subtype_colors


########################################################################################################################
############################################### Import Data from Krug ##################################################

def kr_ImportData():
    """
    Import protein, mRNA, metadata from Krug et al
    """
    # annotation file
    kr_annotation = pd.read_csv('myapp/Data/Krug/kr_group_key.txt', sep='\t')
    kr_annotation = kr_annotation.set_index('Patient')  # set annotation file index as sample name

    # protein file
    kr_df_protein = pd.read_csv('myapp/Data/Krug/kr_data_p.txt', sep='\t')

    # mRNA file
    kr_df_mrna = pd.read_csv('myapp/Data/Krug/kr_data_m.txt', sep='\t')

    # unique gene list
    kr_genes = pd.unique(
        list(kr_df_protein.Gene) + list(kr_df_mrna.Gene))  # combine the genes, this will produce an array
    kr_genes = np.ndarray.tolist(kr_genes)  # convert the array to a list of genes that are not duplicated.

    # transpose dataset
    kr_df_protein = kr_df_protein.rename(columns={'Gene': 'Patient'}).set_index(
        'Patient').T  # rename the column for transposed dataset
    kr_df_mrna = kr_df_mrna.rename(columns={'Gene': 'Patient'}).set_index(
        'Patient').T  # rename the column for transposed dataset

    # subtype information
    kr_subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']  # list of subtypes
    kr_subtype_colors = ['#E31A1C', '#FB9A99', '#1F78B4', '#A6CEE3',
                         '#33A02C']  # list of corresponding colors for subtypes

    kr_df_protein = kr_df_protein.join(kr_annotation.iloc[:, [1]])  # merge subtypes info onto protein dataset
    kr_df_mrna = kr_df_mrna.join(kr_annotation.iloc[:, [1]])  # merge subtypes info onto mRNA dataset

    # sort data by subtype
    kr_df_protein = kr_df_protein.sort_index().sort_values(by='PAM50', ascending=True, kind='mergesort')
    kr_df_mrna = kr_df_mrna.sort_index().sort_values(by='PAM50', ascending=True, kind='mergesort')

    return kr_df_protein, kr_df_mrna, kr_genes, kr_subtypes, kr_subtype_colors


########################################################################################################################
############################################### Import Data from Mertins ###############################################

def me_ImportData():
    """
    Import protein, mRNA, metadata from Mertins et al
    """
    # annotation file
    me_annotation = pd.read_csv('myapp/Data/Mertins/me_group_key.txt', sep='\t')
    me_annotation = me_annotation.set_index('Patient')  # set annotation file index as sample name

    # protein file
    me_df_protein = pd.read_csv('myapp/Data/Mertins/me_data_p.txt', sep='\t')

    # mRNA file
    me_df_mrna = pd.read_csv('myapp/Data/Mertins/me_data_m.txt', sep='\t')


    # unique gene list
    me_genes = pd.unique(
        list(me_df_protein.Gene) + list(me_df_mrna.Gene))  # combine the genes, this will produce an array
    me_genes = np.ndarray.tolist(me_genes)  # convert the array to a list of genes that are not duplicated.

    # transpose dataset
    me_df_protein = me_df_protein.rename(columns={'Gene': 'Patient'}).set_index(
        'Patient').T  # rename the column for transposed dataset
    me_df_mrna = me_df_mrna.rename(columns={'Gene': 'Patient'}).set_index(
        'Patient').T  # rename the column for transposed dataset

    # subtype information
    me_subtypes = ['Basal', 'Her2', 'LumA', 'LumB']  # list of subtypes
    me_subtype_colors = ['#E31A1C', '#FB9A99', '#1F78B4', '#A6CEE3']  # list of corresponding colors for subtypes

    me_df_protein = me_df_protein.join(me_annotation)  # merge subtypes info onto protein dataset
    me_df_mrna = me_df_mrna.join(me_annotation)  # merge subtypes info onto mRNA dataset

    # sort data by subtype
    me_df_protein = me_df_protein.sort_index().sort_values(by='PAM50', ascending=True, kind='mergesort')
    me_df_mrna = me_df_mrna.sort_index().sort_values(by='PAM50', ascending=True, kind='mergesort')

    return me_df_protein, me_df_mrna, me_genes, me_subtypes, me_subtype_colors

def Import_Correlation_Data():
    jo_protein_for_cor = pd.read_csv('myapp/Data/Johansson/jo_data_p.txt',index_col='Gene', sep='\t')
    jo_mrna_for_cor = pd.read_csv('myapp/Data/Johansson/jo_data_m.txt',index_col='Gene', sep='\t')

    kr_protein_for_cor = pd.read_csv('myapp/Data/Krug/kr_data_p.txt',index_col='Gene', sep='\t')
    kr_mrna_for_cor = pd.read_csv('myapp/Data/Krug/kr_data_m.txt',index_col='Gene', sep='\t')

    me_protein_for_cor = pd.read_csv('myapp/Data/Mertins/me_data_p.txt',index_col='Gene', sep='\t')
    me_mrna_for_cor = pd.read_csv('myapp/Data/Mertins/me_data_m.txt',index_col='Gene', sep='\t')

    return jo_mrna_for_cor, kr_mrna_for_cor, me_mrna_for_cor, jo_protein_for_cor, kr_protein_for_cor, me_protein_for_cor