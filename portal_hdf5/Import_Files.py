import pandas as pd


def Import_HDF5():
    JohanssonProteome = '/Users/zhuoheng/Desktop/Vacanti/RawData_March/hdf5/JohanssonProteome.hdf5'
    JohanssonTranscriptome = '/Users/zhuoheng/Desktop/Vacanti/RawData_March/hdf5/JohanssonTranscriptome.hdf5'

    KrugProteome = '/Users/zhuoheng/Desktop/Vacanti/RawData_March/hdf5/KrugProteome.hdf5'
    KrugTranscriptome = '/Users/zhuoheng/Desktop/Vacanti/RawData_March/hdf5/KrugTranscriptome.hdf5'

    MertinsProteome = '/Users/zhuoheng/Desktop/Vacanti/RawData_March/hdf5/MertinsProteome.hdf5'
    MertinsTranscriptome = '/Users/zhuoheng/Desktop/Vacanti/RawData_March/hdf5/MertinsTranscriptome.hdf5'

    return JohanssonProteome, JohanssonTranscriptome, KrugProteome, KrugTranscriptome, MertinsProteome, MertinsTranscriptome


def Import_Static_Correlation_Table():
    jo_mrna_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/jo_mrna_for_cor_ERBB2.txt', sep='\t')
    jo_protein_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/jo_protein_for_cor_ERBB2.txt', sep='\t')

    kr_mrna_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/kr_mrna_for_cor_ERBB2.txt', sep='\t')
    kr_protein_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/kr_protein_for_cor_ERBB2.txt', sep='\t')

    me_mrna_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/me_mrna_for_cor_ERBB2.txt', sep='\t')
    me_protein_ERBB2 = pd.read_csv('portal/Data/StaticCorrelationTable/me_protein_for_cor_ERBB2.txt', sep='\t')

    return jo_mrna_ERBB2, jo_protein_ERBB2, kr_mrna_ERBB2, kr_protein_ERBB2, me_mrna_ERBB2, me_protein_ERBB2