import numpy as np
import Import_data

########################################################################################################################
##################################################### Johansson ########################################################
def jo_Subtype_Avg_SEM_DFs(jo_df_protein, jo_df_mrna, jo_subtypes):
    """
    Calculate the mean and SEM protein and mRNA expression of each gene based on subtypes
    """
    grouped_protein_df = jo_df_protein.groupby(['Pam50']) # group by subtype first
    grouped_mrna_df = jo_df_mrna.groupby(['Pam50'])

    jo_protein_mean_by_subtype = grouped_protein_df.mean() # calculate MEAN of each gene based on subtype
    jo_mrna_mean_by_subtype = grouped_mrna_df.mean()

    jo_protein_sem_by_subtype = grouped_protein_df.sem() # calculate SEM of each gene based on subtype
    jo_mrna_sem_by_subtype = grouped_mrna_df.sem()

    dfs = [jo_protein_mean_by_subtype, jo_mrna_mean_by_subtype, jo_protein_sem_by_subtype, jo_mrna_sem_by_subtype]
    for i in dfs:
        zero_data = np.zeros(len(jo_subtypes))
        i.loc[:, 'blank 1'] = zero_data
        i.loc[:, 'blank 2'] = zero_data
        i.loc[:, 'blank 3'] = zero_data
        i.loc[:, 'blank 4'] = zero_data  # create a 'blank' gene column with zeroes as input

    return jo_protein_mean_by_subtype, jo_mrna_mean_by_subtype, jo_protein_sem_by_subtype, jo_mrna_sem_by_subtype

########################################################################################################################
##################################################### Krug #############################################################

def kr_Subtype_Avg_SEM_DFs(kr_df_protein, kr_df_mrna, kr_subtypes):
    """
    Calculate the mean and SEM protein and mRNA expression of each gene based on subtypes
    """
    grouped_protein_df = kr_df_protein.groupby(['PAM50'])  # group by subtype first
    grouped_mrna_df = kr_df_mrna.groupby(['PAM50'])

    kr_protein_mean_by_subtype = grouped_protein_df.mean()  # calculate MEAN of each gene based on subtype
    kr_mrna_mean_by_subtype = grouped_mrna_df.mean()

    kr_protein_sem_by_subtype = grouped_protein_df.sem()  # calculate SEM of each gene based on subtype
    kr_mrna_sem_by_subtype = grouped_mrna_df.sem()

    dfs = [kr_protein_mean_by_subtype, kr_mrna_mean_by_subtype, kr_protein_sem_by_subtype, kr_mrna_sem_by_subtype]
    for i in dfs:
        zero_data = np.zeros(len(kr_subtypes))
        i.loc[:, 'blank 1'] = zero_data
        i.loc[:, 'blank 2'] = zero_data
        i.loc[:, 'blank 3'] = zero_data
        i.loc[:, 'blank 4'] = zero_data  # create a 'blank' gene column with zeroes as input

    return kr_protein_mean_by_subtype, kr_mrna_mean_by_subtype, kr_protein_sem_by_subtype, kr_mrna_sem_by_subtype

########################################################################################################################
##################################################### Mertins ##########################################################

def me_Subtype_Avg_SEM_DFs(me_df_protein, me_df_mrna, me_subtypes):
    """
    Calculate the mean and SEM protein and mRNA expression of each gene based on subtypes
    """
    grouped_protein_df = me_df_protein.groupby(['PAM50'])  # group by subtype first
    grouped_mrna_df = me_df_mrna.groupby(['PAM50'])

    me_protein_mean_by_subtype = grouped_protein_df.mean()  # calculate MEAN of each gene based on subtype
    me_mrna_mean_by_subtype = grouped_mrna_df.mean()

    me_protein_sem_by_subtype = grouped_protein_df.sem()  # calculate SEM of each gene based on subtype
    me_mrna_sem_by_subtype = grouped_mrna_df.sem()

    dfs = [me_protein_mean_by_subtype, me_mrna_mean_by_subtype, me_protein_sem_by_subtype, me_mrna_sem_by_subtype]
    for i in dfs:
        zero_data = np.zeros(len(me_subtypes))
        i.loc[:, 'blank 1'] = zero_data
        i.loc[:, 'blank 2'] = zero_data
        i.loc[:, 'blank 3'] = zero_data
        i.loc[:, 'blank 4'] = zero_data  # create a 'blank' gene column with zeroes as input

    return me_protein_mean_by_subtype, me_mrna_mean_by_subtype, me_protein_sem_by_subtype, me_mrna_sem_by_subtype

#test
#print(me_Subtype_Avg_SEM_DFs(Import_data.me_ImportData()[0],Import_data.me_ImportData()[1], Import_data.me_ImportData()[3]))







