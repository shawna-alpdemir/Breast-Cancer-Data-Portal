# function to import HDF5 files and static correlation tables

import pandas as pd

def Import_HDF5():
    JohanssonProteome = 'portal/Data/HDF5/JohanssonProteome.hdf5'
    JohanssonTranscriptome = 'portal/Data/HDF5/JohanssonTranscriptome.hdf5'

    KrugProteome = 'portal/Data/HDF5/KrugProteome.hdf5'
    KrugTranscriptome = 'portal/Data/HDF5/KrugTranscriptome.hdf5'

    MertinsProteome = 'portal/Data/HDF5/MertinsProteome.hdf5'
    MertinsTranscriptome = 'portal/Data/HDF5/MertinsTranscriptome.hdf5'

    return JohanssonProteome, JohanssonTranscriptome, KrugProteome, KrugTranscriptome, MertinsProteome, MertinsTranscriptome


def Import_Static_Correlation_Table():
    jo_mrna_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/jo_mrna_for_cor_ERBB2.txt', sep='\t')
    jo_protein_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/jo_protein_for_cor_ERBB2.txt', sep='\t')

    kr_mrna_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/kr_mrna_for_cor_ERBB2.txt', sep='\t')
    kr_protein_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/kr_protein_for_cor_ERBB2.txt', sep='\t')

    me_mrna_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/me_mrna_for_cor_ERBB2.txt', sep='\t')
    me_protein_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/me_protein_for_cor_ERBB2.txt', sep='\t')

    summary_mrna_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/ERBB2_mRNA_summary_table.txt', sep='\t', index_col='Gene')
    summary_protein_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/ERBB2_protein_summary_table.txt', sep='\t', index_col='Gene')

    return jo_mrna_ERBB2, jo_protein_ERBB2, kr_mrna_ERBB2, kr_protein_ERBB2, me_mrna_ERBB2, me_protein_ERBB2, summary_protein_ERBB2, summary_mrna_ERBB2

def Import_Cor_Matrix_HDF5():
    jo_protein_cormat = 'portal/Data/HDF5_Cormat/jo_protein_cormat.hdf5'
    jo_mrna_cormat = 'portal/Data/HDF5_Cormat/jo_mrna_cormat.hdf5'

    kr_protein_cormat = 'portal/Data/HDF5_Cormat/kr_protein_cormat.hdf5'
    kr_mrna_cormat = 'portal/Data/HDF5_Cormat/kr_mrna_cormat.hdf5'

    me_protein_cormat = 'portal/Data/HDF5_Cormat/me_protein_cormat.hdf5'
    me_mrna_cormat = 'portal/Data/HDF5_Cormat/me_mrna_cormat.hdf5'

    annotate = pd.read_csv('portal/Data/Gene name annotation.csv', index_col='Gene')

    return jo_protein_cormat, jo_mrna_cormat, kr_protein_cormat, kr_mrna_cormat, me_protein_cormat, me_mrna_cormat, annotate

def Import_Subtype_DF():
    jo_pro = pd.read_csv('portal/Data/Subtype_avg_sem/jo_pro_subtype.csv', index_col='Gene')
    jo_rna = pd.read_csv('portal/Data/Subtype_avg_sem/jo_rna_subtype.csv', index_col='Gene')

    kr_pro = pd.read_csv('portal/Data/Subtype_avg_sem/kr_pro_subtype.csv', index_col='Gene')
    kr_rna = pd.read_csv('portal/Data/Subtype_avg_sem/kr_rna_subtype.csv', index_col='Gene')

    me_pro = pd.read_csv('portal/Data/Subtype_avg_sem/me_pro_subtype.csv', index_col='Gene')
    me_rna = pd.read_csv('portal/Data/Subtype_avg_sem/me_rna_subtype.csv', index_col='Gene')

    return jo_pro, jo_rna, kr_rna, kr_pro, me_rna, me_pro