import h5py
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
import itertools
from brokenaxes import brokenaxes
import string
import pdb
import matplotlib.pyplot as plt

from Import_Files import Import_HDF5, Import_Subtype_DF

# constant
[JohanssonProteome, JohanssonTranscriptome, KrugProteome, KrugTranscriptome, MertinsProteome, MertinsTranscriptome] = Import_HDF5()
[jo_pro, jo_rna, kr_rna, kr_pro, me_rna, me_pro] = Import_Subtype_DF()
color_list = ['red','blue','green','orange','purple']

# line plot
def Make_Line_Plot(gene_list, dataset, study, yaxis_label, output_filename):

    df_col_name = ['Gene']
    subtype_name = ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']

    if study == 'jo':
        sample_subtype_name = [i + ' ' + str(j) for i in subtype_name for j in range(1, 10)]
        df_col_name = df_col_name + sample_subtype_name
        DF = pd.DataFrame(columns=df_col_name)
    else:
        kr_basal = ['Basal' + ' ' + str(j) for j in range(1, 30)]
        kr_her2 = ['Her2' + ' ' + str(j) for j in range(1, 15)]
        kr_luma = ['LumA' + ' ' + str(j) for j in range(1, 58)]
        kr_lumb = ['LumB' + ' ' + str(j) for j in range(1, 18)]
        kr_norm = ['Norm' + ' ' + str(j) for j in range(1, 6)]
        sample_subtype_name = kr_basal + kr_her2 + kr_luma + kr_lumb + kr_norm
        df_col_name = df_col_name + sample_subtype_name
        DF = pd.DataFrame(columns=df_col_name)

    with h5py.File(dataset, "r") as a:
        for i in gene_list:
            abundance = list(a[i])
            abundance.insert(0, i)
            DF.loc[len(DF.index)] = abundance

    DF.set_index('Gene', inplace=True)
    DF.columns.name = 'Subtype'
    DF_array = pd.DataFrame(DF.stack(), columns=['abundance']).reset_index()

    sns.despine()
    plt.figure(figsize=(6,3))
    plot = sns.lineplot(data=DF_array, x='Subtype', y='abundance', hue='Gene',
                        palette=color_list[:len(gene_list)], dashes=False)
    plt.xlabel('', fontsize=12)
    plt.ylabel(yaxis_label, fontsize=12)
    if study=='jo':
        plt.xticks(['Basal 9','Her2 9', 'LumA 9', 'LumB 9', 'Norm 9'], ['Basal', 'Her2', 'LumA', 'LumB', 'Norm'])
    else:
        plt.xticks(['Basal 29', 'Her2 14', 'LumA 57', 'LumB 17', 'Norm 5'], ['Basal', 'Her2', 'LumA', 'LumB', 'Norm'])
    plt.savefig(f'{output_filename}.pdf',dpi=1500)

def Make_Scat_Plot(gene_list, dataset, study, protein_oder_mrna, output_filename):

    df_col_name = ['Subtype']
    subtype_name = ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']
    colorcode = {'Basal':'#E31A1C', 'Her2': '#FB9A99', 'LumA':'#1F78B4', 'LumB':'#A6CEE3', 'Norm':'#33A02C'}
    if study == 'jo':
        sample_subtype_name = list(itertools.chain.from_iterable(itertools.repeat(x, 9) for x in subtype_name))
        df_col_name = df_col_name + sample_subtype_name
        DF = pd.DataFrame(columns=df_col_name)
    elif study == 'kr':
        sample_subtype_name = ['Basal']*29 + ['Her2']*14 + ['LumA']*57 + ['LumB']*17 + ['Norm']*5
        df_col_name = df_col_name + sample_subtype_name
        DF = pd.DataFrame(columns=df_col_name)
    else:
        sample_subtype_name = ['Basal'] * 18 + ['Her2'] * 12 + ['LumA'] * 23 + ['LumB'] * 24
        df_col_name = df_col_name + sample_subtype_name
        DF = pd.DataFrame(columns=df_col_name)
        colorcode = {'Basal':'#E31A1C', 'Her2': '#FB9A99', 'LumA':'#1F78B4', 'LumB':'#A6CEE3'}

    with h5py.File(dataset, "r") as a:
        for i in gene_list:
            abundance = list(a[i])
            abundance.insert(0, i)
            DF.loc[len(DF.index)] = abundance

    DF.set_index('Subtype', inplace=True)
    DF_array = DF.transpose()

    # plot except 2B 2C 2D
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x=DF_array[gene_list[0]], y=DF_array[gene_list[1]], c=DF_array.index.map(colorcode)) # scatter plot with color mapping
    plt.xlabel(f'{gene_list[0]} {protein_oder_mrna} z-score', fontsize=12) # x axis label
    plt.ylabel(f'{gene_list[1]} {protein_oder_mrna} z-score', fontsize=12) # y axis label

    #  fig 2A 2B 2C 2D specifically # ylims=((-1.5, 2), (4.5, 5.5)),
    # fig = plt.figure(figsize=(5,4))
    # bax = brokenaxes(xlims=((-0.5, 1), (6.25, 6.5)),
    #                  hspace=.05)
    #
    # #scatterplots
    # bax.scatter(x=DF_array[gene_list[0]], y=DF_array[gene_list[1]],c=DF_array.index.map(colorcode))  # scatter plot with color mapping
    # bax.set_xlabel(f'{gene_list[0]} {protein_oder_mrna} z-score', fontsize=12) # x axis label
    # bax.set_ylabel(f'{gene_list[1]} {protein_oder_mrna} z-score', fontsize=12) # y axis label

    # call the scipy function for pearson correlation
    r, p = scipy.stats.pearsonr(x=DF_array[gene_list[0]], y = DF_array[gene_list[1]])
    # Plot regression line
    ax.plot(DF_array[gene_list[0]], p + r * DF_array[gene_list[0]], color="k", lw=1)

    # annotate the pearson correlation coefficient text to 2 decimal places
    plt.text(.05, .8, 'r={:.2f}'.format(r))
    plt.text(.05, .75, 'p={:.3}'.format(p))
    plt.savefig(f'{output_filename}.pdf',dpi=1500)

def Make_Subtype_Plot(gene_list, dataset, study, protein_oder_mrna, output_filename):
    # johansson and krug: study == 20, mertins: study == 16
    avg_DF = pd.DataFrame()
    sem_DF = pd.DataFrame()
    for i in gene_list:
         with h5py.File(dataset, "r") as a:
            subtype_list = np.array(a.get('subtypes'))
            subtype_list = [x.decode('utf-8') for x in subtype_list]

            gene = np.array(a.get(i))

            gene_DF = pd.DataFrame(data={'subtype': subtype_list, 'values': gene})

            avg_DF = avg_DF.append(gene_DF.groupby(by='subtype').mean())
            sem_DF = sem_DF.append(gene_DF.groupby(by='subtype').sem())

            avg_DF.fillna(0,inplace=True)
            sem_DF.fillna(0,inplace=True)

    avg_DF['sem'] = sem_DF['values']

    avg_DF['alphabet'] = list(string.ascii_lowercase[:study])

    plt.errorbar(data=avg_DF,x='alphabet',y='values',yerr='sem', linestyle='', fmt='o', markersize=5, capsize=5)
    plt.xticks('')
    plt.xlabel('')
    plt.ylabel(f'{protein_oder_mrna} z-score')
    plt.savefig(f'{output_filename}.pdf', dpi=1500)

# -----------------------------------x
# Figure 1, S1
# -----------------------------------
# 1A
# Make_Line_Plot(gene_list=['NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS7', 'NDUFS8'],
#                dataset=JohanssonProteome,
#                study='jo',
#                yaxis_label='protein z-score',
#                output_filename='Johansson_complexI_pro')
# S1A
# Make_Line_Plot(gene_list=['NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS7', 'NDUFS8'],
#                dataset=KrugProteome,
#                study='kr',
#                yaxis_label='protein z-score',
#                output_filename='Krug_S1A')
# 1B
# Make_Line_Plot(gene_list=['NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS7', 'NDUFS8'],
#                dataset=JohanssonTranscriptome,
#                study='jo',
#                yaxis_label='mRNA z-score',
#                output_filename='Johansson_complexI_mRNA')
# S1B
# Make_Line_Plot(gene_list=['NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS7', 'NDUFS8'],
#                dataset=KrugTranscriptome,
#                study='kr',
#                yaxis_label='mRNA z-score',
#                output_filename='Krug_S1B')
# 1C
# Make_Scat_Plot(gene_list=['PDHA1','PDHB'],
#                dataset=JohanssonProteome,
#                study='jo',
#                protein_oder_mrna='protein',
#                output_filename='Johansson_PDH_complex')
# S1C
# Make_Scat_Plot(gene_list=['PDHA1','PDHB'],
#                dataset=KrugProteome,
#                study='kr',
#                protein_oder_mrna='protein',
#                output_filename='Krug_S1C')
# 1D
# Make_Scat_Plot(gene_list=['SDHA','SDHB'],
#                dataset=JohanssonProteome,
#                study='jo',
#                protein_oder_mrna='protein',
#                output_filename='Johansson_SDH_complex')
# S1D
# Make_Scat_Plot(gene_list=['SDHA','SDHB'],
#                dataset=KrugProteome,
#                study='kr',
#                protein_oder_mrna='protein',
#                output_filename='Krug_S1D')
# 1E
# Make_Scat_Plot(gene_list=['CYC1','UQCRH'],
#                dataset=JohanssonProteome,
#                study='jo',
#                protein_oder_mrna='protein',
#                output_filename='Johansson_cytC_red_complex')
# S1E
# Make_Scat_Plot(gene_list=['CYC1','UQCRH'],
#                dataset=KrugProteome,
#                study='kr',
#                protein_oder_mrna='protein',
#                output_filename='Krug_S1E')
# 1F
# Make_Line_Plot(gene_list=['COX5A','COX5B','COX7A2'],
#                dataset=JohanssonProteome,
#                study='jo',
#                yaxis_label= 'protein z-score',
#                output_filename='Johansson_cytC_ox_complex')
# S1F
Make_Line_Plot(gene_list=['COX5A','COX5B','COX7A2'],
               dataset=KrugProteome,
               study='kr',
               yaxis_label= 'protein z-score',
               output_filename='Krug_S1F')

# -----------------------------------
# Figure 2, S2
# -----------------------------------
# 2A
# Make_Scat_Plot(gene_list=['SHMT2','SLC25A32'],
#                dataset=JohanssonProteome,
#                study='jo',
#                protein_oder_mrna='protein',
#                output_filename='Johansson_fig2A')
# 2B
# Make_Scat_Plot(gene_list=['SHMT2','MTHFR'],
#                dataset=JohanssonProteome,
#                study='jo',
#                protein_oder_mrna='protein',
#                output_filename='Johansson_fig2B')
# 2C
# Make_Scat_Plot(gene_list=['PHGDH','SHMT2'],
#                dataset=JohanssonProteome,
#                study='jo',
#                protein_oder_mrna='protein',
#                output_filename='Johansson_fig2C')
# 2D
# Make_Scat_Plot(gene_list=['PHGDH','SHMT1'],
#                dataset=JohanssonProteome,
#                study='jo',
#                protein_oder_mrna='protein',
#                output_filename='Johansson_fig2D')
# 2E
# Make_Subtype_Plot(gene_list=['SHMT2','SLC25A32','MTHFR','PHGDH'],
#                   dataset=JohanssonProteome,
#                   study=20,
#                   protein_oder_mrna='protein',
#                   output_filename='Johansson_protein_2E1')
# Make_Subtype_Plot(gene_list=['SHMT2','SLC25A32','MTHFR','PHGDH'],
#                   dataset=JohanssonTranscriptome,
#                   study=20,
#                   protein_oder_mrna='mRNA',
#                   output_filename='Johansson_mRNA_2E2')
# S2A
# Make_Scat_Plot(gene_list=['PHGDH','SHMT2'],
#                dataset=KrugProteome,
#                study='kr',
#                protein_oder_mrna='protein',
#                output_filename='Krug_figS2A1')
# Make_Scat_Plot(gene_list=['PHGDH','SHMT1'],
#                dataset=KrugProteome,
#                study='kr',
#                protein_oder_mrna='protein',
#                output_filename='Krug_figS2A2')
# S2B
# Make_Scat_Plot(gene_list=['PHGDH','SHMT2'],
#                dataset=MertinsProteome,
#                study='me',
#                protein_oder_mrna='protein',
#                output_filename='Mertins_figS2B1')
# Make_Scat_Plot(gene_list=['PHGDH','SHMT1'],
#                dataset=MertinsProteome,
#                study='me',
#                protein_oder_mrna='protein',
#                output_filename='Mertins_figS2B2')
# S2C
# Make_Subtype_Plot(gene_list=['SHMT2','SLC25A32','MTHFR','PHGDH'],
#                   dataset=KrugProteome,
#                   study=20,
#                   protein_oder_mrna='protein',
#                   output_filename='Krug_protein_S2C')
# Make_Subtype_Plot(gene_list=['SHMT2','SLC25A32','MTHFR','PHGDH'],
#                   dataset=KrugTranscriptome,
#                   study=20,
#                   protein_oder_mrna='mRNA',
#                   output_filename='Krug_mRNA_S2C')
#
# Make_Subtype_Plot(gene_list=['SHMT2','SLC25A32','MTHFR','PHGDH'],
#                   dataset=MertinsProteome,
#                   study=16,
#                   protein_oder_mrna='protein',
#                   output_filename='Mertins_protein_S2C')
# Make_Subtype_Plot(gene_list=['SHMT2','SLC25A32','MTHFR','PHGDH'],
#                   dataset=MertinsTranscriptome,
#                   study=16,
#                   protein_oder_mrna='mRNA',
#                   output_filename='Mertins_mRNA_S2C')