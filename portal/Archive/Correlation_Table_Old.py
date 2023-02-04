# Archived Feb 4, 2023
# Old method: use Proteome and Transcriptome HDF5 file to create correlation matrix from scratch

# def Get_Protein_Correlation_Table(gene_name):
#
#     with h5py.File(JohanssonProteome, "r") as a:
#
#         try:
#             # input gene
#             gene = np.array(a.get(gene_name))
#
#             # total genes
#             total_gene_list = np.array(a.get('genes'))
#             total_gene_list = [x.decode('utf-8') for x in total_gene_list]
#             total_gene_list.remove(gene_name)  # list without input gene
#
#             array = []
#             for i in total_gene_list:
#                 array.append(np.array(a.get(i)))
#             # gene_df = pd.DataFrame(data=gene, index=['ERBB2']
#             rest_df = pd.DataFrame(data=array, index=total_gene_list).dropna()
#             gene_df = pd.Series(data=gene).transpose()
#
#             Correlation_matrix = rest_df.corrwith(gene_df, axis=1)  # perform correlation
#             Correlation_array = np.asarray(Correlation_matrix)
#             Correlation_array = np.round(Correlation_array, 4)  # round r value
#
#             TestStat = (Correlation_matrix * (45 - 2) ** 0.5) / (1 - Correlation_matrix ** 2) ** 0.5  # calculate t score
#             Sig = sc.t.sf(abs(TestStat),df=45 - 2) * 2  # using the survival function to compute one sided p-value, then doubled it to make it two-side
#
#             Gene_name_array = np.asarray(rest_df.index, dtype=str)  # get the gene name
#             Correlation_table_df = pd.DataFrame(data=np.column_stack((Gene_name_array, Correlation_array, Sig)),columns=['Gene', 'r', 'p'], index=None)
#             Correlation_table_df[['r', 'p']] = Correlation_table_df[['r', 'p']].apply(pd.to_numeric)
#             Correlation_table_df = Correlation_table_df.sort_values('p',ascending=True)  # sort the gene by smallest p value
#             jo_Correlation_table_df = Correlation_table_df.head(100)
#
#         except ValueError:
#             jo_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
#                                                     "r": ["NA"],
#                                                     "p": ["NA"]}, index=None)
#
#     with h5py.File(KrugProteome, "r") as a:
#
#         try:
#             # input gene
#             gene = np.array(a.get(gene_name))
#
#             # total genes
#             total_gene_list = np.array(a.get('genes'))
#             total_gene_list = [x.decode('utf-8') for x in total_gene_list]
#             total_gene_list.remove(gene_name) # list without input gene
#
#             array = []
#             for i in total_gene_list:
#                 array.append(np.array(a.get(i)))
#
#             rest_df = pd.DataFrame(data=array, index=total_gene_list).dropna()
#             gene_df = pd.Series(data=gene).transpose()
#
#             Correlation_matrix = rest_df.corrwith(gene_df, axis=1)  # perform correlation
#             Correlation_array = np.asarray(Correlation_matrix)
#             Correlation_array = np.round(Correlation_array, 4)  # round r value
#
#             TestStat = (Correlation_matrix * (122 - 2) ** 0.5) / (1 - Correlation_matrix ** 2) ** 0.5  # calculate t score
#             Sig = sc.t.sf(abs(TestStat), df=122 - 2) * 2  # using the survival function to compute one sided p-value, then doubled it to make it two-side
#
#             Gene_name_array = np.asarray(rest_df.index, dtype=str)  # get the gene name
#             Correlation_table_df = pd.DataFrame(data=np.column_stack((Gene_name_array, Correlation_array, Sig)),columns=['Gene', 'r', 'p'], index=None)
#             Correlation_table_df[['r', 'p']] = Correlation_table_df[['r', 'p']].apply(pd.to_numeric)
#             Correlation_table_df = Correlation_table_df.sort_values('p',ascending=True)  # sort the gene by smallest p value
#             kr_Correlation_table_df = Correlation_table_df.head(100)
#
#         except ValueError:
#             kr_Correlation_table_df = pd.DataFrame({"Gene" : [f"No data for {gene_name}"],
#                                             "r": ["NA"],
#                                             "p": ["NA"]}, index=None)
#
#     with h5py.File(MertinsProteome, "r") as a:
#
#         try:
#             # input gene
#             gene = np.array(a.get(gene_name))
#
#              # total genes
#             total_gene_list = np.array(a.get('genes'))
#             total_gene_list = [x.decode('utf-8') for x in total_gene_list]
#             total_gene_list.remove(gene_name)  # list without input gene
#
#             array = []
#             for i in total_gene_list:
#                 array.append(np.array(a.get(i)))
#
#             rest_df = pd.DataFrame(data=array, index=total_gene_list).dropna()
#             gene_df = pd.Series(data=gene).transpose()
#
#             Correlation_matrix = rest_df.corrwith(gene_df, axis=1)  # perform correlation
#             Correlation_array = np.asarray(Correlation_matrix)
#             Correlation_array = np.round(Correlation_array, 4)  # round r value
#
#             TestStat = (Correlation_matrix * (77 - 2) ** 0.5) / (1 - Correlation_matrix ** 2) ** 0.5  # calculate t score
#             Sig = sc.t.sf(abs(TestStat), df=77 - 2) * 2  # using the survival function to compute one sided p-value, then doubled it to make it two-side
#
#             Gene_name_array = np.asarray(rest_df.index, dtype=str)  # get the gene name
#             Correlation_table_df = pd.DataFrame(data=np.column_stack((Gene_name_array, Correlation_array, Sig)),columns=['Gene', 'r', 'p'], index=None)
#             Correlation_table_df[['r', 'p']] = Correlation_table_df[['r', 'p']].apply(pd.to_numeric)
#             Correlation_table_df = Correlation_table_df.sort_values('p',ascending=True)  # sort the gene by smallest p value
#             me_Correlation_table_df = Correlation_table_df.head(100)
#
#         except ValueError:
#             me_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
#                                                     "r": ["NA"],
#                                                     "p": ["NA"]}, index=None)
#     return jo_Correlation_table_df, kr_Correlation_table_df, me_Correlation_table_df
#
# def Get_mRNA_Correlation_Table(gene_name):
#
#     with h5py.File(JohanssonTranscriptome, "r") as a:
#
#         try:
#             # input gene
#             gene = np.array(a.get(gene_name))
#
#             # total genes
#             total_gene_list = np.array(a.get('genes'))
#             total_gene_list = [x.decode('utf-8') for x in total_gene_list]
#             total_gene_list.remove(gene_name)  # list without input gene
#
#             array = []
#             for i in total_gene_list:
#                 array.append(np.array(a.get(i)))
#             # gene_df = pd.DataFrame(data=gene, index=['ERBB2']
#             rest_df = pd.DataFrame(data=array, index=total_gene_list).dropna()
#             gene_df = pd.Series(data=gene).transpose()
#
#             Correlation_matrix = rest_df.corrwith(gene_df, axis=1)  # perform correlation
#             Correlation_array = np.asarray(Correlation_matrix)
#             Correlation_array = np.round(Correlation_array, 4)  # round r value
#
#             TestStat = (Correlation_matrix * (45 - 2) ** 0.5) / (1 - Correlation_matrix ** 2) ** 0.5  # calculate t score
#             Sig = sc.t.sf(abs(TestStat),df=45 - 2) * 2  # using the survival function to compute one sided p-value, then doubled it to make it two-side
#
#             Gene_name_array = np.asarray(rest_df.index, dtype=str)  # get the gene name
#             Correlation_table_df = pd.DataFrame(data=np.column_stack((Gene_name_array, Correlation_array, Sig)),columns=['Gene', 'r', 'p'], index=None)
#             Correlation_table_df[['r', 'p']] = Correlation_table_df[['r', 'p']].apply(pd.to_numeric)
#             Correlation_table_df = Correlation_table_df.sort_values('p',ascending=True)  # sort the gene by smallest p value
#             jo_m_Correlation_table_df = Correlation_table_df.head(100)
#
#         except ValueError:
#             jo_m_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
#                                                     "r": ["NA"],
#                                                     "p": ["NA"]}, index=None)
#
#     with h5py.File(KrugTranscriptome, "r") as a:
#
#         try:
#             # input gene
#             gene = np.array(a.get(gene_name))
#
#             # total genes
#             total_gene_list = np.array(a.get('genes'))
#             total_gene_list = [x.decode('utf-8') for x in total_gene_list]
#             total_gene_list.remove(gene_name) # list without input gene
#
#             array = []
#             for i in total_gene_list:
#                 array.append(np.array(a.get(i)))
#
#             rest_df = pd.DataFrame(data=array, index=total_gene_list).dropna()
#             gene_df = pd.Series(data=gene).transpose()
#
#             Correlation_matrix = rest_df.corrwith(gene_df, axis=1)  # perform correlation
#             Correlation_array = np.asarray(Correlation_matrix)
#             Correlation_array = np.round(Correlation_array, 4)  # round r value
#
#             TestStat = (Correlation_matrix * (122 - 2) ** 0.5) / (1 - Correlation_matrix ** 2) ** 0.5  # calculate t score
#             Sig = sc.t.sf(abs(TestStat), df=122 - 2) * 2  # using the survival function to compute one sided p-value, then doubled it to make it two-side
#
#             Gene_name_array = np.asarray(rest_df.index, dtype=str)  # get the gene name
#             Correlation_table_df = pd.DataFrame(data=np.column_stack((Gene_name_array, Correlation_array, Sig)),columns=['Gene', 'r', 'p'], index=None)
#             Correlation_table_df[['r', 'p']] = Correlation_table_df[['r', 'p']].apply(pd.to_numeric)
#             Correlation_table_df = Correlation_table_df.sort_values('p',ascending=True)  # sort the gene by smallest p value
#             kr_m_Correlation_table_df = Correlation_table_df.head(100)
#
#         except ValueError:
#             kr_m_Correlation_table_df = pd.DataFrame({"Gene" : [f"No data for {gene_name}"],
#                                             "r": ["NA"],
#                                             "p": ["NA"]}, index=None)
#
#     with h5py.File(MertinsTranscriptome, "r") as a:
#
#         try:
#             # input gene
#             gene = np.array(a.get(gene_name))
#
#              # total genes
#             total_gene_list = np.array(a.get('genes'))
#             total_gene_list = [x.decode('utf-8') for x in total_gene_list]
#             total_gene_list.remove(gene_name)  # list without input gene
#
#             array = []
#             for i in total_gene_list:
#                 array.append(np.array(a.get(i)))
#
#             rest_df = pd.DataFrame(data=array, index=total_gene_list).dropna()
#             gene_df = pd.Series(data=gene).transpose()
#
#             Correlation_matrix = rest_df.corrwith(gene_df, axis=1)  # perform correlation
#             Correlation_array = np.asarray(Correlation_matrix)
#             Correlation_array = np.round(Correlation_array, 4)  # round r value
#
#             TestStat = (Correlation_matrix * (77 - 2) ** 0.5) / (1 - Correlation_matrix ** 2) ** 0.5  # calculate t score
#             Sig = sc.t.sf(abs(TestStat), df=77 - 2) * 2  # using the survival function to compute one sided p-value, then doubled it to make it two-side
#
#             Gene_name_array = np.asarray(rest_df.index, dtype=str)  # get the gene name
#             Correlation_table_df = pd.DataFrame(data=np.column_stack((Gene_name_array, Correlation_array, Sig)),columns=['Gene', 'r', 'p'], index=None)
#             Correlation_table_df[['r', 'p']] = Correlation_table_df[['r', 'p']].apply(pd.to_numeric)
#             Correlation_table_df = Correlation_table_df.sort_values('p',ascending=True)  # sort the gene by smallest p value
#             me_m_Correlation_table_df = Correlation_table_df.head(100)
#
#         except ValueError:
#             me_m_Correlation_table_df = pd.DataFrame({"Gene": [f"No data for {gene_name}"],
#                                                     "r": ["NA"],
#                                                     "p": ["NA"]}, index=None)
#     return jo_m_Correlation_table_df, kr_m_Correlation_table_df, me_m_Correlation_table_df
