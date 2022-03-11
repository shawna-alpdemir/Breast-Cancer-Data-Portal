import pandas as pd
import scipy.stats as sc
import numpy as np
import time

jo_mrna_for_cor = pd.read_csv('/Users/zhuoheng/PycharmProjects/Breast Cancer Data Portal/portal/Data/Johansson/jo_data_m.txt',
                              index_col='Gene', sep='\t')
jo_protein_for_cor = pd.read_csv('/Users/zhuoheng/PycharmProjects/Breast Cancer Data Portal/portal/Data/Johansson/jo_data_p.txt',
                                 index_col='Gene', sep='\t')

kr_mrna_for_cor = pd.read_csv('/Users/zhuoheng/PycharmProjects/Breast Cancer Data Portal/portal/Data/Krug/kr_data_m.txt',index_col='Gene', sep='\t')
kr_protein_for_cor = pd.read_csv('/Users/zhuoheng/PycharmProjects/Breast Cancer Data Portal/portal/Data/Krug/kr_data_p.txt',index_col='Gene', sep='\t')

me_mrna_for_cor = pd.read_csv('/Users/zhuoheng/PycharmProjects/Breast Cancer Data Portal/portal/Data/Mertins/me_data_m.txt',index_col='Gene', sep='\t')
me_protein_for_cor = pd.read_csv('/Users/zhuoheng/PycharmProjects/Breast Cancer Data Portal/portal/Data/Mertins/me_data_p.txt',index_col='Gene', sep='\t')

list = [jo_mrna_for_cor,jo_protein_for_cor,kr_mrna_for_cor,kr_protein_for_cor,me_mrna_for_cor,me_protein_for_cor]
name_list = ['jo_mrna_for_cor','jo_protein_for_cor','kr_mrna_for_cor','kr_protein_for_cor','me_mrna_for_cor','me_protein_for_cor']


def GeneCorrelationTable(Dataset, Gene_name, Output_name):
    """ generate static table for any genes """
    # start_time = time.time()

    Dataset_without_input = Dataset.drop(Gene_name)
    # drop the gene in the textbox, so the input gene can be correlated with the rest of the genes.

    Correlation_matrix = Dataset_without_input.transpose().corrwith(
        Dataset.transpose().loc[:, Gene_name])  # perform correlation
    Correlation_array = np.asarray(Correlation_matrix)
    Correlation_array = np.round(Correlation_array, 4)  # round r value

    TestStat = (Correlation_matrix * (45 - 2) ** 0.5) / (1 - Correlation_matrix ** 2) ** 0.5  # calculate t score
    Sig = sc.t.sf(abs(TestStat),df=45-2)*2  # using the survival function to compute one sided p-value, then doubled it to make it two-side

    Gene_name_array = np.asarray(Dataset_without_input.index, dtype=str)  # get the gene name
    Correlation_table_df = pd.DataFrame(data=np.column_stack((Gene_name_array, Correlation_array, Sig)),
                                        columns=['Gene', 'r', 'p'])
    Correlation_table_df = Correlation_table_df.sort_values('p', ascending=True)  # sort the gene by smallest p value
    Correlation_table_df = Correlation_table_df.head(100)

    Correlation_table_df.to_csv(f'/Users/zhuoheng/PycharmProjects/Breast Cancer Data Portal/portal/Data/StaticCorrelationTable/{Output_name}.txt', sep='\t', index=False)

    # end_time = time.time()
    # time_elapsed = end_time - start_time
    # print(time_elapsed)

for i,j in zip(list,name_list):
    GeneCorrelationTable(i, 'ERBB2', j)