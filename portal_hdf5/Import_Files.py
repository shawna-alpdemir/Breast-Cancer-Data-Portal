# function to import HDF5 files and static correlation tables

import pandas as pd

def Import_HDF5():
    JohanssonProteome = 'portal_hdf5/Data/HDF5/JohanssonProteome.hdf5'
    JohanssonTranscriptome = 'portal_hdf5/Data/HDF5/JohanssonTranscriptome.hdf5'

    KrugProteome = 'portal_hdf5/Data/HDF5/KrugProteome.hdf5'
    KrugTranscriptome = 'portal_hdf5/Data/HDF5/KrugTranscriptome.hdf5'

    MertinsProteome = 'portal_hdf5/Data/HDF5/MertinsProteome.hdf5'
    MertinsTranscriptome = 'portal_hdf5/Data/HDF5/MertinsTranscriptome.hdf5'

    return JohanssonProteome, JohanssonTranscriptome, KrugProteome, KrugTranscriptome, MertinsProteome, MertinsTranscriptome


def Import_Static_Correlation_Table():
    jo_mrna_ERBB2 = pd.read_csv('portal_hdf5/Data/StaticCorrelationTable/jo_mrna_for_cor_ERBB2.txt', sep='\t')
    jo_protein_ERBB2 = pd.read_csv('portal_hdf5/Data/StaticCorrelationTable/jo_protein_for_cor_ERBB2.txt', sep='\t')

    kr_mrna_ERBB2 = pd.read_csv('portal_hdf5/Data/StaticCorrelationTable/kr_mrna_for_cor_ERBB2.txt', sep='\t')
    kr_protein_ERBB2 = pd.read_csv('portal_hdf5/Data/StaticCorrelationTable/kr_protein_for_cor_ERBB2.txt', sep='\t')

    me_mrna_ERBB2 = pd.read_csv('portal_hdf5/Data/StaticCorrelationTable/me_mrna_for_cor_ERBB2.txt', sep='\t')
    me_protein_ERBB2 = pd.read_csv('portal_hdf5/Data/StaticCorrelationTable/me_protein_for_cor_ERBB2.txt', sep='\t')

    return jo_mrna_ERBB2, jo_protein_ERBB2, kr_mrna_ERBB2, kr_protein_ERBB2, me_mrna_ERBB2, me_protein_ERBB2