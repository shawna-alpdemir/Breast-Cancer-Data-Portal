# Archived Feb 4, 2023. Heatmap may be a feature in the future

import h5py
import numpy as np
import pandas as pd
import pdb as pdb
from bokeh.models import BasicTicker, ColorBar, LinearColorMapper
from bokeh.plotting import figure, show
from Correlation_Table import Get_Protein_Correlation_Table, Get_mRNA_Correlation_Table
from Import_Files import Import_HDF5
from Plot_All_Styling import StylePlots
# import abundance data
[JohanssonProteome, JohanssonTranscriptome, KrugProteome, KrugTranscriptome, MertinsProteome, MertinsTranscriptome] = Import_HDF5()

# set up color map, color using Sunset diverging colour scheme from https://personal.sron.nl/~pault/. this palette is related to RdYlBu scheme
colors = ['#364B9A', '#4A7BB7', '#6EA6CD', '#98CAE1', '#C2E4EF', '#EAECCC', '#FEDA8B', '#FDB366', '#F67E4B', '#DD3D2D', '#A50026']
mapper = LinearColorMapper(palette=colors, low=-5, high=5)

# function call, and get the top genes in summary table.
def Get_Protein_Heatmaps():

    pro_cor_tables = pd.read_csv('portal_hdf5/Data/StaticCorrelationTable/erbb2_pro_summed_table.txt',sep='\t',index_col='Gene')
    pro_sum_gene_list = list(pro_cor_tables.head(30).index)

    # construct column names
    df_col_name = ['Gene']
    subtype_name = ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']

    jo_sample_subtype_name = [i + ' ' + str(j) for i in subtype_name for j in range(1, 10)]

    kr_basal = ['Basal'+' '+str(j) for j in range(1,30)]
    kr_her2 = ['Her2' + ' ' + str(j) for j in range(1, 15)]
    kr_luma = ['LumA' + ' ' + str(j) for j in range(1, 58)]
    kr_lumb = ['LumB' + ' ' + str(j) for j in range(1, 18)]
    kr_norm = ['Norm' + ' ' + str(j) for j in range(1, 6)]
    kr_sample_subtype_name = kr_basal+kr_her2+kr_luma+kr_lumb+kr_norm

    me_basal = ['Basal'+' '+str(j) for j in range(1,19)]
    me_her2 = ['Her2' + ' ' + str(j) for j in range(1, 13)]
    me_luma = ['LumA' + ' ' + str(j) for j in range(1, 24)]
    me_lumb = ['LumB' + ' ' + str(j) for j in range(1, 25)]
    me_sample_subtype_name = me_basal+me_her2+me_luma+me_lumb

    jo_df_col_name = df_col_name + jo_sample_subtype_name
    kr_df_col_name = df_col_name + kr_sample_subtype_name
    me_df_col_name = df_col_name + me_sample_subtype_name

    jo_pro_hm_df = pd.DataFrame(columns=jo_df_col_name)
    kr_pro_hm_df = pd.DataFrame(columns=kr_df_col_name)
    me_pro_hm_df = pd.DataFrame(columns=me_df_col_name)

    with h5py.File(JohanssonProteome, "r") as a:

    # if Nate wants tumor ID as col names
    # jo_df_col_name = ['Gene']
    # tumor_list = np.array(a.get('tumors'))
    # tumor_list = [x.decode('utf-8') for x in tumor_list]
    # jo_df_col_name.extend(tumor_list)
    # jo_pro_hm_df = pd.DataFrame(columns=jo_df_col_name)

        for i in pro_sum_gene_list:
            try:
                abundance = list(a[i])
                abundance.insert(0, i)
                jo_pro_hm_df.loc[len(jo_pro_hm_df.index)] = abundance
            except KeyError:
                pass

    with h5py.File(KrugProteome, "r") as a:
        for i in pro_sum_gene_list:
            try:
                abundance = list(a[i])
                abundance.insert(0, i)
                kr_pro_hm_df.loc[len(kr_pro_hm_df.index)] = abundance
            except KeyError:
                pass

    with h5py.File(MertinsProteome, "r") as a:
        for i in pro_sum_gene_list:
            try:
                abundance = list(a[i])
                abundance.insert(0, i)
                me_pro_hm_df.loc[len(me_pro_hm_df.index)] = abundance
            except KeyError:
                pass

    # set up index and columns for each study
    jo_pro_hm_df.set_index('Gene', inplace=True)
    jo_pro_hm_df.columns.name = 'Subtype'
    jo_pro_genes = list(jo_pro_hm_df.index)

    kr_pro_hm_df.set_index('Gene', inplace=True)
    kr_pro_hm_df.columns.name = 'Subtype'
    kr_pro_genes = list(kr_pro_hm_df.index)

    me_pro_hm_df.set_index('Gene', inplace=True)
    me_pro_hm_df.columns.name = 'Subtype'
    me_pro_genes = list(me_pro_hm_df.index)

    # reshape to 1D array or rates with a gene and subtype for each row.
    jo_pro_array = pd.DataFrame(jo_pro_hm_df.stack(), columns=['abundance']).reset_index()
    kr_pro_array = pd.DataFrame(kr_pro_hm_df.stack(), columns=['abundance']).reset_index()
    me_pro_array = pd.DataFrame(me_pro_hm_df.stack(), columns=['abundance']).reset_index()

    # set up heatmaps
    jo_pro_hm = figure(title="Johansson heatmap for top 30 correlated genes".format(jo_sample_subtype_name[-1], jo_sample_subtype_name[0]),
           x_range=jo_sample_subtype_name, y_range=list(reversed(jo_pro_genes)),
           x_axis_location="above", width=1100, height=400,
           tooltips=[('gene', '@Gene'), ('subtype', '@Subtype'), ('abundance', '@abundance')])

    kr_pro_hm = figure(title="Krug heatmap for top 30 correlated genes".format(kr_sample_subtype_name[-1], kr_sample_subtype_name[0]),
           x_range=kr_sample_subtype_name, y_range=list(reversed(kr_pro_genes)),
           x_axis_location="above", width=1100, height=400,
           tooltips=[('gene', '@Gene'), ('subtype', '@Subtype'), ('abundance', '@abundance')])

    me_pro_hm = figure(title="Krug heatmap for top 30 correlated genes".format(me_sample_subtype_name[-1], me_sample_subtype_name[0]),
           x_range=me_sample_subtype_name, y_range=list(reversed(me_pro_genes)),
           x_axis_location="above", width=1100, height=400,
           tooltips=[('gene', '@Gene'), ('subtype', '@Subtype'), ('abundance', '@abundance')])

    # set up heatmap figure sources
    jo_pro_hm.rect(x="Subtype", y="Gene", width=1, height=1, source=jo_pro_array,
                      fill_color={'field': 'abundance', 'transform': mapper}, line_color=None)
    kr_pro_hm.rect(x="Subtype", y="Gene", width=1, height=1, source=kr_pro_array,
                      fill_color={'field': 'abundance', 'transform': mapper}, line_color=None)
    me_pro_hm.rect(x="Subtype", y="Gene", width=1, height=1, source=me_pro_array,
                      fill_color={'field': 'abundance', 'transform': mapper}, line_color=None)

    # styling the heatmap
    jo_pro_hm = StylePlots(jo_pro_hm, PlotID='Heatmap')
    kr_pro_hm = StylePlots(kr_pro_hm, PlotID='Heatmap')
    me_pro_hm = StylePlots(me_pro_hm, PlotID='Heatmap')

    return jo_pro_hm, kr_pro_hm, me_pro_hm


def Get_mRNA_Heatmaps():

    mRNA_cor_tables = pd.read_csv('portal_hdf5/Data/StaticCorrelationTable/erbb2_mrna_summed_table.txt', sep='\t', index_col='Gene')
    mRNA_sum_gene_list = mRNA_cor_tables.head(30).index

    # construct column names
    df_col_name = ['Gene']
    subtype_name = ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']

    jo_sample_subtype_name = [i + ' ' + str(j) for i in subtype_name for j in range(1, 10)]

    kr_basal = ['Basal'+' '+str(j) for j in range(1,30)]
    kr_her2 = ['Her2' + ' ' + str(j) for j in range(1, 15)]
    kr_luma = ['LumA' + ' ' + str(j) for j in range(1, 58)]
    kr_lumb = ['LumB' + ' ' + str(j) for j in range(1, 18)]
    kr_norm = ['Norm' + ' ' + str(j) for j in range(1, 6)]
    kr_sample_subtype_name = kr_basal+kr_her2+kr_luma+kr_lumb+kr_norm

    me_basal = ['Basal'+' '+str(j) for j in range(1,19)]
    me_her2 = ['Her2' + ' ' + str(j) for j in range(1, 13)]
    me_luma = ['LumA' + ' ' + str(j) for j in range(1, 24)]
    me_lumb = ['LumB' + ' ' + str(j) for j in range(1, 25)]
    me_sample_subtype_name = me_basal+me_her2+me_luma+me_lumb

    jo_df_col_name = df_col_name + jo_sample_subtype_name
    kr_df_col_name = df_col_name + kr_sample_subtype_name
    me_df_col_name = df_col_name + me_sample_subtype_name

    jo_mRNA_hm_df = pd.DataFrame(columns=jo_df_col_name)
    kr_mRNA_hm_df = pd.DataFrame(columns=kr_df_col_name)
    me_mRNA_hm_df = pd.DataFrame(columns=me_df_col_name)


    with h5py.File(JohanssonTranscriptome, "r") as a:

        # if Nate wants tumor ID as col names
        # jo_df_col_name = ['Gene']
        # tumor_list = np.array(a.get('tumors'))
        # tumor_list = [x.decode('utf-8') for x in tumor_list]
        # jo_df_col_name.extend(tumor_list)
        # jo_pro_hm_df = pd.DataFrame(columns=jo_df_col_name)

        for i in mRNA_sum_gene_list:
            try:
                abundance = list(a[i])
                abundance.insert(0, i)
                jo_mRNA_hm_df.loc[len(jo_mRNA_hm_df.index)] = abundance
            except KeyError:
                pass

    with h5py.File(KrugTranscriptome, "r") as a:
        for i in mRNA_sum_gene_list:
            try:
                abundance = list(a[i])
                abundance.insert(0, i)
                kr_mRNA_hm_df.loc[len(kr_mRNA_hm_df.index)] = abundance
            except KeyError:
                pass

    with h5py.File(MertinsTranscriptome, "r") as a:
        for i in mRNA_sum_gene_list:
            try:
                abundance = list(a[i])
                abundance.insert(0, i)
                me_mRNA_hm_df.loc[len(me_mRNA_hm_df.index)] = abundance
            except KeyError:
                pass

    # set up index and columns for each study
    jo_mRNA_hm_df.set_index('Gene', inplace=True)
    jo_mRNA_hm_df.columns.name = 'Subtype'
    jo_mRNA_genes = list(jo_mRNA_hm_df.index)

    kr_mRNA_hm_df.set_index('Gene', inplace=True)
    kr_mRNA_hm_df.columns.name = 'Subtype'
    kr_mRNA_genes = list(kr_mRNA_hm_df.index)

    me_mRNA_hm_df.set_index('Gene', inplace=True)
    me_mRNA_hm_df.columns.name = 'Subtype'
    me_mRNA_genes = list(me_mRNA_hm_df.index)

    # reshape to 1D array or rates with a gene and subtype for each row.
    jo_mRNA_array = pd.DataFrame(jo_mRNA_hm_df.stack(), columns=['abundance']).reset_index()
    kr_mRNA_array = pd.DataFrame(kr_mRNA_hm_df.stack(), columns=['abundance']).reset_index()
    me_mRNA_array = pd.DataFrame(me_mRNA_hm_df.stack(), columns=['abundance']).reset_index()

    # set up heatmaps
    jo_mRNA_hm = figure(title="Johansson heatmap for top 30 correlated genes".format(jo_sample_subtype_name[-1], jo_sample_subtype_name[0]),
                       x_range=jo_sample_subtype_name, y_range=list(reversed(jo_mRNA_genes)),
                       x_axis_location="above", width=1100, height=400,
                        tooltips=[('gene', '@Gene'), ('subtype', '@Subtype'), ('abundance', '@abundance')])

    kr_mRNA_hm = figure(title="Krug heatmap for top 30 correlated genes".format(kr_sample_subtype_name[-1], kr_sample_subtype_name[0]),
                       x_range=kr_sample_subtype_name, y_range=list(reversed(kr_mRNA_genes)),
                       x_axis_location="above", width=1100, height=400,
                        tooltips=[('gene', '@Gene'), ('subtype', '@Subtype'), ('abundance', '@abundance')])

    me_mRNA_hm = figure(title="Mertins heatmap for top 30 correlated genes".format(me_sample_subtype_name[-1], me_sample_subtype_name[0]),
                       x_range=me_sample_subtype_name, y_range=list(reversed(me_mRNA_genes)),
                       x_axis_location="above", width=1100, height=400,
                        tooltips=[('gene', '@Gene'), ('subtype', '@Subtype'), ('abundance', '@abundance')])

    # set up heatmap figure sources
    jo_mRNA_hm.rect(x="Subtype", y="Gene", width=1, height=1, source=jo_mRNA_array,
                      fill_color={'field': 'abundance', 'transform': mapper}, line_color=None)
    kr_mRNA_hm.rect(x="Subtype", y="Gene", width=1, height=1, source=kr_mRNA_array,
                      fill_color={'field': 'abundance', 'transform': mapper}, line_color=None)
    me_mRNA_hm.rect(x="Subtype", y="Gene", width=1, height=1, source=me_mRNA_array,
                      fill_color={'field': 'abundance', 'transform': mapper}, line_color=None)

    # styling the heatmap
    jo_mRNA_hm = StylePlots(jo_mRNA_hm, PlotID='Heatmap')
    kr_mRNA_hm = StylePlots(kr_mRNA_hm, PlotID='Heatmap')
    me_mRNA_hm = StylePlots(me_mRNA_hm, PlotID='Heatmap')

    return jo_mRNA_hm, kr_mRNA_hm, me_mRNA_hm