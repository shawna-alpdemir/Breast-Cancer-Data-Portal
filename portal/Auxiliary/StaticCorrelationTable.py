import pandas as pd
import scipy.stats as sc
import numpy as np

jo_mrna_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/jo_mrna_z.csv', index_col='Gene')
jo_protein_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/jo_pro_z.csv', index_col='Gene')

kr_mrna_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/kr_rna_z.csv',index_col='Gene')
kr_protein_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/kr_pro_z.csv', index_col='Gene')

me_mrna_for_cor = pd .read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/me_rna_z.csv', index_col='Gene')
me_protein_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/CleanData_April/me_pro_z.csv', index_col='Gene')

annotate = pd.read_csv('/Users/zhuoheng/PycharmProjects/Breast Cancer Data Portal/portal_hdf5/Data/Gene name annotation.csv', index_col='Gene')

list = [jo_mrna_for_cor,jo_protein_for_cor,kr_mrna_for_cor,kr_protein_for_cor,me_mrna_for_cor,me_protein_for_cor]
name_list = ['jo_mrna_for_cor_ERBB2','jo_protein_for_cor_ERBB2','kr_mrna_for_cor_ERBB2','kr_protein_for_cor_ERBB2','me_mrna_for_cor_ERBB2','me_protein_for_cor_ERBB2']


# only for proteome data, filter out genes that have more than 30% of NA
# perc = 30.0 # percentage of allowing NA
# min_count =  int(((100-perc)/100)*jo_protein_for_cor.shape[1] + 1)
# mod_df = jo_protein_for_cor.dropna(axis=0,
#                     thresh=min_count)
# print(jo_protein_for_cor.shape)
# print(mod_df.shape)

def GeneCorrelationTable(Dataset, Gene_name, Output_name):
    """ generate static table for any genes """
    # start_time = time.time()

    Dataset_without_input = Dataset.drop(Gene_name)
    # drop the gene in the textbox, so the input gene can be correlated with the rest of the genes.

    Correlation_matrix = Dataset_without_input.transpose().corrwith(
        Dataset.transpose().loc[:, Gene_name])  # perform correlation
    Correlation_array = np.asarray(Correlation_matrix)
    Correlation_array = np.round(Correlation_array, 4)  # round r value

    TestStat = (Correlation_matrix * (77 - 2) ** 0.5) / (1 - Correlation_matrix ** 2) ** 0.5  # calculate t score
    Sig = sc.t.sf(abs(TestStat),df=77-2)*2  # using the survival function to compute one sided p-value, then doubled it to make it two-side

    Gene_name_array = np.asarray(Dataset_without_input.index, dtype=str)  # get the gene name
    Correlation_table_df = pd.DataFrame(data=np.column_stack((Gene_name_array, Correlation_array, Sig)),
                                        columns=['Gene', 'r', 'p'])
    Correlation_table_df = Correlation_table_df.join(annotate, on='Gene')
    #Correlation_table_df[['r','p']] = Correlation_table_df[['r','p']].apply(pd.to_numeric)
    Correlation_table_df['r'] = pd.to_numeric(Correlation_table_df['r'], errors='coerce')
    Correlation_table_df['p'] = pd.to_numeric(Correlation_table_df['p'], errors='coerce')
    Correlation_table_df = Correlation_table_df.sort_values('p', ascending=True)  # sort the gene by smallest p value
    Correlation_table_df = Correlation_table_df.head(1000)
    Correlation_table_df.to_csv(f'/Users/zhuoheng/PycharmProjects/Breast Cancer Data Portal/portal_hdf5/Data/StaticCorrelationTable/{Output_name}.txt', sep='\t', index=False)

    # end_time = time.time()
    # time_elapsed = end_time - start_time
    # print(time_elapsed)



#GeneCorrelationTable(jo_mrna_for_cor,'ERBB2','jo_mrna_for_cor_ERBB2')
for i,j in zip([me_mrna_for_cor,me_protein_for_cor],['me_mrna_for_cor_ERBB2','me_protein_for_cor_ERBB2']):
   GeneCorrelationTable(i, 'ERBB2', j)
