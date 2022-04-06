# function that creates subtype plot dataframes for means and SEMs

import h5py
import numpy as np
import pandas as pd
from pdb import set_trace
from Import_Files import Import_HDF5

# function call
[JohanssonProteome, JohanssonTranscriptome, KrugProteome, KrugTranscriptome, MertinsProteome, MertinsTranscriptome] = Import_HDF5()

########################################################################################################################
##################################################### Johansson ########################################################
def Check_NaN_in_DF(df):
    """
    Check if the DF contains NaN, if so, swap NaN with 0
    """
    is_nan_exist = df.isnull().any()
    if is_nan_exist[0]: # return true then fill with 0, false then pass
        df.fillna(0, inplace=True)
    else:
        pass

def Johansson_Subtype_Avg_SEM_DFs(gene):
    """
    Calculate the mean and SEM protein and mRNA expression of each gene based on subtypes
    """
    with h5py.File(JohanssonProteome, "r") as a, h5py.File(JohanssonTranscriptome, "r") as b:
        subtype_list = np.array(a.get('subtypes'))
        subtype_list = [x.decode('utf-8') for x in subtype_list]

        gene_protein = np.array(a.get(gene))
        gene_mRNA = np.array(b.get(gene))

        gene_protein_DF = pd.DataFrame(data={'subtype': subtype_list, 'values': gene_protein})
        gene_mRNA_DF = pd.DataFrame(data={'subtype': subtype_list, 'values': gene_mRNA})

        jo_avg_protein_DF = gene_protein_DF.groupby(by='subtype').mean()
        jo_sem_protein_DF = gene_protein_DF.groupby(by='subtype').sem()

        jo_avg_mRNA_DF = gene_mRNA_DF.groupby(by='subtype').mean()
        jo_sem_mRNA_DF = gene_mRNA_DF.groupby(by='subtype').sem()

        Check_NaN_in_DF(jo_avg_protein_DF)
        Check_NaN_in_DF(jo_sem_protein_DF)
        Check_NaN_in_DF(jo_avg_mRNA_DF)
        Check_NaN_in_DF(jo_sem_mRNA_DF)

    return jo_avg_protein_DF, jo_sem_protein_DF, jo_avg_mRNA_DF, jo_sem_mRNA_DF


########################################################################################################################
##################################################### Krug #############################################################

def Krug_Subtype_Avg_SEM_DFs(gene):
    """
    Calculate the mean and SEM protein and mRNA expression of each gene based on subtypes
    """
    with h5py.File(KrugProteome, "r") as a, h5py.File(KrugTranscriptome, "r") as b:
        subtype_list = np.array(a.get('subtypes'))
        subtype_list = [x.decode('utf-8') for x in subtype_list]

        gene_protein = np.array(a.get(gene))
        gene_mRNA = np.array(b.get(gene))

        gene_protein_DF = pd.DataFrame(data={'subtype': subtype_list, 'values': gene_protein})
        gene_mRNA_DF = pd.DataFrame(data={'subtype': subtype_list, 'values': gene_mRNA})

        kr_avg_protein_DF = gene_protein_DF.groupby(by='subtype').mean()
        kr_sem_protein_DF = gene_protein_DF.groupby(by='subtype').sem()

        kr_avg_mRNA_DF = gene_mRNA_DF.groupby(by='subtype').mean()
        kr_sem_mRNA_DF = gene_mRNA_DF.groupby(by='subtype').sem()

        Check_NaN_in_DF(kr_avg_protein_DF)
        Check_NaN_in_DF(kr_sem_protein_DF)
        Check_NaN_in_DF(kr_avg_mRNA_DF)
        Check_NaN_in_DF(kr_sem_mRNA_DF)

    return kr_avg_protein_DF, kr_sem_protein_DF, kr_avg_mRNA_DF, kr_sem_mRNA_DF


########################################################################################################################
##################################################### Mertins ##########################################################

def Mertins_Subtype_Avg_SEM_DFs(gene):
    """
    Calculate the mean and SEM protein and mRNA expression of each gene based on subtypes
    """
    with h5py.File(MertinsProteome, "r") as a, h5py.File(MertinsTranscriptome, "r") as b:
        subtype_list = np.array(a.get('subtypes'))
        subtype_list = [x.decode('utf-8') for x in subtype_list]

        gene_protein = np.array(a.get(gene))
        gene_mRNA = np.array(b.get(gene))

        gene_protein_DF = pd.DataFrame(data={'subtype': subtype_list, 'values': gene_protein})
        gene_mRNA_DF = pd.DataFrame(data={'subtype': subtype_list, 'values': gene_mRNA})

        me_avg_protein_DF = gene_protein_DF.groupby(by='subtype').mean()
        me_sem_protein_DF = gene_protein_DF.groupby(by='subtype').sem()

        me_avg_mRNA_DF = gene_mRNA_DF.groupby(by='subtype').mean()
        me_sem_mRNA_DF = gene_mRNA_DF.groupby(by='subtype').sem()

        Check_NaN_in_DF(me_avg_protein_DF)
        Check_NaN_in_DF(me_sem_protein_DF)
        Check_NaN_in_DF(me_avg_mRNA_DF)
        Check_NaN_in_DF(me_sem_mRNA_DF)

    return me_avg_protein_DF, me_sem_protein_DF, me_avg_mRNA_DF, me_sem_mRNA_DF