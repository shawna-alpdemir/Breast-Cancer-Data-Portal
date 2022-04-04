import pandas as pd
import scipy.stats as sc
import numpy as np
import h5py
import time

#jo_mrna_for_cor = pd.read_csv('/portal/Data/Johansson/jo_data_m.txt',index_col='Gene', sep='\t')
jo_protein_for_cor = pd.read_csv('/Users/zhuoheng/Desktop/Vacanti/RawData_March/ConvertedDateCleanData/jo_data_p.txt',index_col='Gene', sep='\t')

#kr_mrna_for_cor = pd.read_csv('/portal/Data/Krug/kr_data_m.txt', index_col='Gene', sep='\t')
#kr_protein_for_cor = pd.read_csv('/portal/Data/Krug/kr_data_p.txt', index_col='Gene', sep='\t')

#me_mrna_for_cor = pd.read_csv('/portal/Data/Mertins/me_data_m.txt', index_col='Gene', sep='\t')
#me_protein_for_cor = pd.read_csv('/portal/Data/Mertins/me_data_p.txt', index_col='Gene', sep='\t')

#list = [jo_mrna_for_cor,jo_protein_for_cor,kr_mrna_for_cor,kr_protein_for_cor,me_mrna_for_cor,me_protein_for_cor]
#name_list = ['jo_mrna_for_cor_ERBB2','jo_protein_for_cor_ERBB2','kr_mrna_for_cor_ERBB2','kr_protein_for_cor_ERBB2','me_mrna_for_cor_ERBB2','me_protein_for_cor_ERBB2']


JohanssonProteinCorrelation = '/Users/zhuoheng/Desktop/Vacanti/RawData_March/hdf5/JohanssonProteinCorrelation.hdf5'

def GeneCorrelationTable(Dataset, Gene_name):
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
    Correlation_table_df[['r','p']] = Correlation_table_df[['r','p']].apply(pd.to_numeric)
    Correlation_table_df = Correlation_table_df.sort_values('p', ascending=True)  # sort the gene by smallest p value
    Correlation_table_df = Correlation_table_df.head(100)

    #Correlation_table_array = Correlation_table_df.to_numpy()
    df_list = []
    for column in Correlation_table_df.columns:
        data = Correlation_table_df[column].tolist()
        df_list.append(data)

    return df_list
    #Correlation_table_df.to_csv(f'/Users/zhuoheng/PycharmProjects/Breast Cancer Data Portal/portal/Data/StaticCorrelationTable/{Output_name}.txt', sep='\t', index=False)

    # end_time = time.time()
    # time_elapsed = end_time - start_time
    # print(time_elapsed)

start = time.time()
print("start")
with h5py.File(JohanssonProteinCorrelation, "w") as hdf5:
    #for i in jo_protein_for_cor.index:
    df_list = GeneCorrelationTable(jo_protein_for_cor, 'ERBB2')
    hdf5.create_dataset('ERBB2', data=np.array(df_list, dtype='S'))
end = time.time()
elapse = end-start
print(elapse)

with h5py.File(JohanssonProteinCorrelation, "r") as hdf5:
    erbb2_data = np.array(hdf5.get('ERBB2'))
    erbb2_data = [x.decode('utf-8') for x in erbb2_data[0]]
print(erbb2_data)