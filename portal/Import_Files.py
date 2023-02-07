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

def Import_Coef_Matrix_HDF5():
    jo_protein_coef = 'portal/Data/HDF5_Cormat/jo_pro_coef.hdf5'
    jo_mrna_coef = 'portal/Data/HDF5_Cormat/jo_mrna_coef.hdf5'

    kr_protein_cormat = 'portal/Data/HDF5_Cormat/kr_pro_coef.hdf5'
    kr_mrna_coef = 'portal/Data/HDF5_Cormat/kr_mrna_coef.hdf5'

    me_protein_coef = 'portal/Data/HDF5_Cormat/me_pro_coef.hdf5'
    me_mrna_coef = 'portal/Data/HDF5_Cormat/me_mrna_coef.hdf5'

    annotate = pd.read_csv('portal/Data/Gene name annotation.csv', index_col='Gene')

    return jo_protein_coef, jo_mrna_coef, kr_protein_cormat, kr_mrna_coef, me_protein_coef, me_mrna_coef, annotate

def Import_Pval_Matrix_HDF5():
    jo_protein_pval = 'portal/Data/HDF5_Cormat/jo_pro_pval.hdf5'
    jo_mrna_pval = 'portal/Data/HDF5_Cormat/jo_mrna_pval.hdf5'

    kr_protein_pval = 'portal/Data/HDF5_Cormat/kr_pro_pval.hdf5'
    kr_mrna_pval = 'portal/Data/HDF5_Cormat/kr_mrna_pval.hdf5'

    me_protein_pval = 'portal/Data/HDF5_Cormat/me_pro_pval.hdf5'
    me_mrna_pval = 'portal/Data/HDF5_Cormat/me_mrna_pval.hdf5'

    return jo_protein_pval, jo_mrna_pval, kr_protein_pval, kr_mrna_pval, me_protein_pval, me_mrna_pval

def Import_Subtype_DF():
    jo_pro = pd.read_csv('portal/Data/Subtype_avg_sem/jo_pro_subtype.csv', index_col='Gene')
    jo_rna = pd.read_csv('portal/Data/Subtype_avg_sem/jo_rna_subtype.csv', index_col='Gene')

    kr_pro = pd.read_csv('portal/Data/Subtype_avg_sem/kr_pro_subtype.csv', index_col='Gene')
    kr_rna = pd.read_csv('portal/Data/Subtype_avg_sem/kr_rna_subtype.csv', index_col='Gene')

    me_pro = pd.read_csv('portal/Data/Subtype_avg_sem/me_pro_subtype.csv', index_col='Gene')
    me_rna = pd.read_csv('portal/Data/Subtype_avg_sem/me_rna_subtype.csv', index_col='Gene')

    return jo_pro, jo_rna, kr_rna, kr_pro, me_rna, me_pro