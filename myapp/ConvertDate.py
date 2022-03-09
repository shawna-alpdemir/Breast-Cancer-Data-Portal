import pandas as pd
import numpy as np

def ConvertDates(DF, name):
    """
    Revert the gene symbols that were mistakenly converted to date by Excel
    """
    df_conv = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/GeneDateConversions.txt', sep='\t')
    df_conv.set_index('Date', inplace=True, drop=False)
    Dates = np.asarray(df_conv['Date'])
    Genes = np.asarray(DF['Gene'])
    for date in Dates:
        if date in Genes:
            DateIndex = np.where(Genes == date)[0]
            NewSymbol = df_conv.loc[date, 'Current']
            Genes[DateIndex] = NewSymbol
    DF['Gene'] = Genes

    DF.to_csv(f'/Users/zhuoheng/Desktop/Vacanti/RawData_March/ConvertedDateCleanData/{name}.txt', sep='\t', index=False)
    return DF

me_rna = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData/me_rna_dropna_ordered_z.txt', sep='\t')
me_protein = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData/me_protein_dropna_ordered_dropdup_z.txt', sep='\t')

kr_rna = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData/kr_rna_dropna_z.txt', sep='\t')
kr_protein = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData/kr_protein_dropna_ordered_dropdup_z.txt', sep='\t')

jo_rna = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData/jo_mrna_dropna_z_rm_nonanno.txt', sep='\t')
jo_protein = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData/jo_protein_dropna_z.txt', sep='\t')


ConvertDates(me_rna, 'me_data_m')
ConvertDates(me_protein, 'me_data_p')

ConvertDates(kr_rna, 'kr_data_m')
ConvertDates(kr_protein, 'kr_data_p')

ConvertDates(jo_rna, 'jo_data_m')
ConvertDates(jo_protein, 'jo_data_p')


