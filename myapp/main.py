######################################### Row 0: Import all necessary packages #########################################
# import python packages
from os.path import join, dirname
import numpy as np
import pandas as pd
import scipy.stats as sc

# import bokeh packages
from bokeh.io import curdoc
from bokeh.layouts import column, row, layout, Spacer
from bokeh.models import FactorRange, AutocompleteInput, Panel, Tabs, ColumnDataSource, DataTable, TableColumn, Button, \
    CustomJS
from bokeh.models.widgets.markups import Div
from bokeh.plotting import figure
from bokeh.transform import factor_cmap

# import necessary functions
from Import_data import jo_ImportData, kr_ImportData, me_ImportData, Import_Correlation_Data
from Correlation_LinePlot_Source import jo_Cor1_Source, kr_Cor1_Source, me_Cor1_Source
from Correlation_CirPlot_Source import jo_Cor2_Source, kr_Cor2_Source, me_Cor2_Source
from Style_Plot import StylePlots
from Subtype_Average_DF import jo_Subtype_Avg_SEM_DFs, kr_Subtype_Avg_SEM_DFs, me_Subtype_Avg_SEM_DFs
from Subtype_Average_Plot_Source import jo_Subtype_Avg_Plot_Source, kr_Subtype_Avg_Plot_Source, me_Subtype_Avg_Plot_Source, reverseTuple
from Subtype_Plot import jo_Subtype_Plot, kr_Subtype_Plot, me_Subtype_Plot

######################################### Row 0: Import all necessary functions ########################################
# From Import_data.py
[jo_df_protein, jo_df_mrna, jo_genes, jo_subtypes, jo_subtype_colors] = jo_ImportData()
[kr_df_protein, kr_df_mrna, kr_genes, kr_subtypes, kr_subtype_colors] = kr_ImportData()
[me_df_protein, me_df_mrna, me_genes, me_subtypes, me_subtype_colors] = me_ImportData()
[jo_mrna_for_cor, kr_mrna_for_cor, me_mrna_for_cor, jo_protein_for_cor, kr_protein_for_cor, me_protein_for_cor] = Import_Correlation_Data()

# From Correlation_linePlot_Source.py
jo_LineSourceList = jo_Cor1_Source(jo_df_protein, jo_df_mrna, INITIAL_GENE=['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4'])
kr_LineSourceList = kr_Cor1_Source(kr_df_protein, kr_df_mrna, INITIAL_GENE=['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4'])
me_LineSourceList = me_Cor1_Source(me_df_protein, me_df_mrna, INITIAL_GENE=['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4'])

# From Correlation_CirPlot_Source.py
jo_cor2_source = jo_Cor2_Source(jo_df_protein, jo_df_mrna, INITIAL_GENE=['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4'])
kr_cor2_source = kr_Cor2_Source(kr_df_protein, kr_df_mrna, INITIAL_GENE=['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4'])
me_cor2_source = me_Cor2_Source(me_df_protein, me_df_mrna, INITIAL_GENE=['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4'])

# From Subtype_Average_DF.py
[jo_protein_mean_by_subtype, jo_mrna_mean_by_subtype, jo_protein_sem_by_subtype,jo_mrna_sem_by_subtype] = jo_Subtype_Avg_SEM_DFs(jo_df_protein, jo_df_mrna, jo_subtypes)
[kr_protein_mean_by_subtype, kr_mrna_mean_by_subtype, kr_protein_sem_by_subtype,kr_mrna_sem_by_subtype] = kr_Subtype_Avg_SEM_DFs(kr_df_protein, kr_df_mrna, kr_subtypes)
[me_protein_mean_by_subtype, me_mrna_mean_by_subtype, me_protein_sem_by_subtype,me_mrna_sem_by_subtype] = me_Subtype_Avg_SEM_DFs(me_df_protein, me_df_mrna, me_subtypes)

# From Subtype_Average_Plot_Source.py
[jo_source_dict_subtype_protein, jo_source_dict_subtype_mrna] = jo_Subtype_Avg_Plot_Source(jo_protein_mean_by_subtype,jo_mrna_mean_by_subtype,jo_protein_sem_by_subtype,jo_mrna_sem_by_subtype,jo_subtypes,gene_for_subtypes_plot=['ESR1', 'PGR', 'ERBB2','MKI67'])
[kr_source_dict_subtype_protein, kr_source_dict_subtype_mrna] = kr_Subtype_Avg_Plot_Source(kr_protein_mean_by_subtype,kr_mrna_mean_by_subtype,kr_protein_sem_by_subtype,kr_mrna_sem_by_subtype,kr_subtypes,gene_for_subtypes_plot=['ESR1', 'PGR', 'ERBB2','MKI67'])
[me_source_dict_subtype_protein, me_source_dict_subtype_mrna] = me_Subtype_Avg_Plot_Source(me_protein_mean_by_subtype,me_mrna_mean_by_subtype,me_protein_sem_by_subtype,me_mrna_sem_by_subtype,me_subtypes,gene_for_subtypes_plot=['ESR1', 'PGR', 'ERBB2','MKI67'])

# From Subtype_Plot.py
[jo_subtype_plot_protein, jo_subtype_plot_mrna] = jo_Subtype_Plot(jo_source_dict_subtype_protein,jo_source_dict_subtype_mrna,jo_subtype_colors,jo_subtypes)
[kr_subtype_plot_protein, kr_subtype_plot_mrna] = kr_Subtype_Plot(kr_source_dict_subtype_protein,kr_source_dict_subtype_mrna,kr_subtype_colors,kr_subtypes)
[me_subtype_plot_protein, me_subtype_plot_mrna] = me_Subtype_Plot(me_source_dict_subtype_protein,me_source_dict_subtype_mrna,me_subtype_colors,me_subtypes)

# Constant
TICKER = {}
TICKER_INDEX = None
TICKER_GENELIST = np.ndarray.tolist(pd.unique(jo_genes + kr_genes + me_genes + ['blank 1','blank 2','blank 3','blank 4']))
INITIAL_GENE = ['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4']
GENE_FOR_SUBTYPES_PLOT = ['ESR1', 'PGR', 'ERBB2', 'MKI67']
GENE_COLORS = ['red', 'blue', 'green', 'orange']
GENE_NUMBER = ['Gene 1', 'Gene 2', 'Gene 3', 'Gene 4']
Gene = 'ERBB2'

####################################################### CALLBACK #######################################################
def nix(val, lst):
    """function that prevents user input the same gene in other tickers"""
    return [x for x in lst if x!= val]


def ticker1_change(attrname, old, new):
    """ change function that triggers the update"""
    print(f"ticker 1 changes to {new}")
    global TICKER_INDEX
    TICKER_INDEX = 0
    TICKER[1].completions = nix(new, TICKER_GENELIST)
    TICKER[2].completions = nix(new, TICKER_GENELIST)
    TICKER[3].completions = nix(new, TICKER_GENELIST)


def ticker2_change(attrname, old, new):
    print(f"ticker 2 changes to {new}")
    global TICKER_INDEX
    TICKER_INDEX = 1
    TICKER[0].completions = nix(new, TICKER_GENELIST)
    TICKER[2].completions = nix(new, TICKER_GENELIST)
    TICKER[3].completions = nix(new, TICKER_GENELIST)


def ticker3_change(attrname, old, new):
    print(f"ticker 3 changes to {new}")
    global TICKER_INDEX
    TICKER_INDEX = 2
    TICKER[0].completions = nix(new, TICKER_GENELIST)
    TICKER[1].completions = nix(new, TICKER_GENELIST)
    TICKER[3].completions = nix(new, TICKER_GENELIST)

def ticker4_change(attrname, old, new):
    print(f"ticker 4 changes to {new}")
    global TICKER_INDEX
    TICKER_INDEX = 3
    TICKER[0].completions = nix(new, TICKER_GENELIST)
    TICKER[1].completions = nix(new, TICKER_GENELIST)
    TICKER[2].completions = nix(new, TICKER_GENELIST)


TICKER_FUNCTION_LIST = [ticker1_change, ticker2_change, ticker3_change, ticker4_change]

def correlation_update(attrname, old, new):
    gene = new
    protein_data = []
    mRNA_data = []
    cor_y = []
    cor_x = []
    mRNA_prot_xList = ['x1', 'x2', 'x3', 'x4']
    mRNA_prot_yList = ['y1', 'y2', 'y3', 'y4']

    for i, j in zip([jo_df_protein, kr_df_protein, me_df_protein], [jo_df_mrna, kr_df_mrna, me_df_mrna]):
        # complex subunit line plots & scatter plot y axis - protein
        if gene in list(i.columns):  # list of gene is now the first row
            protein_data.append(i.loc[:, gene])
            cor_y.append(i.loc[:, gene])
        else:
            protein_data.append(i.loc[:, f'blank {TICKER_INDEX+1}'])
            cor_y.append(i.loc[:, f'blank {TICKER_INDEX+1}'])

        # complex subunit line plots & scatter plot x axis - mRNA
        if gene in list(j.columns):
            mRNA_data.append(j.loc[:, gene])
            cor_x.append(j.loc[:, gene])
        else:
            mRNA_data.append(j.loc[:, f'blank {TICKER_INDEX+1}'])
            cor_x.append(j.loc[:, f'blank {TICKER_INDEX+1}'])

    jo_LineSourceList[TICKER_INDEX].data = dict(subtype=tuple(zip(jo_df_protein['Pam50'], jo_df_protein.index)), protein_data=protein_data[0],
                                                mRNA_data=mRNA_data[0], gene=np.repeat(gene, len(jo_df_protein['Pam50'])))
    kr_LineSourceList[TICKER_INDEX].data = dict(subtype=tuple(zip(kr_df_protein['PAM50'], kr_df_protein.index)), protein_data=protein_data[1],
                                                mRNA_data=mRNA_data[1], gene=np.repeat(gene, len(kr_df_protein['PAM50'])))
    me_LineSourceList[TICKER_INDEX].data = dict(subtype=tuple(zip(me_df_protein['PAM50'], me_df_protein.index)), protein_data=protein_data[2],
                                                mRNA_data=mRNA_data[2], gene=np.repeat(gene, len(me_df_protein['PAM50'])))

    # use ticker to determine which x y is changing
    x_to_change = mRNA_prot_xList[TICKER_INDEX]
    y_to_change = mRNA_prot_yList[TICKER_INDEX]

    jo_cor2_source.data[x_to_change] = cor_x[0]
    jo_cor2_source.data[y_to_change] = cor_y[0]
    kr_cor2_source.data[x_to_change] = cor_x[1]
    kr_cor2_source.data[y_to_change] = cor_y[1]
    me_cor2_source.data[x_to_change] = cor_x[2]
    me_cor2_source.data[y_to_change] = cor_y[2]

    # change line plot legend and scatter plot



def jo_subtype_update(attrname, old, new):
    """call back to update subtype plot for Johansson"""
    gene = new

    x_axis_subtype_protein = []
    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    x_axis_subtype_mrna = []
    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    GENE_FOR_SUBTYPES_PLOT[TICKER_INDEX] = gene # replace the GENE_FOR_SUBTYPES_PLOT list with the new gene
    for gene_to_plot in GENE_FOR_SUBTYPES_PLOT:
        for subtype in jo_subtypes:
            # protein
            x_axis_subtype_protein.append((subtype, gene_to_plot)) # get the gene name for the x axis
            try:
                y_axis_subtype_protein.append((jo_protein_mean_by_subtype.loc[subtype, gene_to_plot]))
                upper_bar_protein.append(jo_protein_mean_by_subtype.loc[subtype, gene_to_plot] + jo_protein_sem_by_subtype.loc[subtype, gene_to_plot])
                lower_bar_protein.append(jo_protein_mean_by_subtype.loc[subtype, gene_to_plot] - jo_protein_sem_by_subtype.loc[subtype, gene_to_plot])
            except KeyError:
                y_axis_subtype_protein.append(0)
                upper_bar_protein.append(0)
                lower_bar_protein.append(0)

            # mRNA
            x_axis_subtype_mrna.append((subtype, gene_to_plot))
            try:
                y_axis_subtype_mrna.append((jo_mrna_mean_by_subtype.loc[subtype, gene_to_plot]))
                upper_bar_mrna.append(jo_mrna_mean_by_subtype.loc[subtype, gene_to_plot] + jo_mrna_sem_by_subtype.loc[subtype, gene_to_plot])
                lower_bar_mrna.append(jo_mrna_mean_by_subtype.loc[subtype, gene_to_plot] - jo_mrna_sem_by_subtype.loc[subtype, gene_to_plot])
            except KeyError:
                y_axis_subtype_mrna.append(0)
                upper_bar_mrna.append(0)
                lower_bar_mrna.append(0)

    # update the plot data
    # protein
    jo_subtype_plot_protein.x_range.factors = reverseTuple(x_axis_subtype_protein)
    jo_source_dict_subtype_protein.data = dict(x=reverseTuple(x_axis_subtype_protein), y=y_axis_subtype_protein,
                                               upper=upper_bar_protein, lower=lower_bar_protein)
    # mRNA
    jo_subtype_plot_mrna.x_range.factors = reverseTuple(x_axis_subtype_mrna)
    jo_source_dict_subtype_mrna.data = dict(x=reverseTuple(x_axis_subtype_mrna), y=y_axis_subtype_mrna,
                                            upper=upper_bar_mrna, lower=lower_bar_mrna)


def kr_subtype_update(attrname, old, new):
    """call back to update subtype plot for Krug"""
    gene = new

    x_axis_subtype_protein = []
    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    x_axis_subtype_mrna = []
    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    GENE_FOR_SUBTYPES_PLOT[TICKER_INDEX] = gene # replace the GENE_FOR_SUBTYPES_PLOT list with the new gene
    for gene_to_plot in GENE_FOR_SUBTYPES_PLOT:
        for subtype in kr_subtypes:
            # protein
            x_axis_subtype_protein.append((subtype, gene_to_plot)) # get the gene name for the x axis
            try:
                y_axis_subtype_protein.append((kr_protein_mean_by_subtype.loc[subtype, gene_to_plot]))
                upper_bar_protein.append(kr_protein_mean_by_subtype.loc[subtype, gene_to_plot] + kr_protein_sem_by_subtype.loc[subtype, gene_to_plot])
                lower_bar_protein.append(kr_protein_mean_by_subtype.loc[subtype, gene_to_plot] - kr_protein_sem_by_subtype.loc[subtype, gene_to_plot])
            except KeyError:
                y_axis_subtype_protein.append(0)
                upper_bar_protein.append(0)
                lower_bar_protein.append(0)

            # mRNA
            x_axis_subtype_mrna.append((subtype, gene_to_plot))
            try:
                y_axis_subtype_mrna.append((kr_mrna_mean_by_subtype.loc[subtype, gene_to_plot]))
                upper_bar_mrna.append(kr_mrna_mean_by_subtype.loc[subtype, gene_to_plot] + kr_mrna_sem_by_subtype.loc[subtype, gene_to_plot])
                lower_bar_mrna.append(kr_mrna_mean_by_subtype.loc[subtype, gene_to_plot] - kr_mrna_sem_by_subtype.loc[subtype, gene_to_plot])
            except KeyError:
                y_axis_subtype_mrna.append(0)
                upper_bar_mrna.append(0)
                lower_bar_mrna.append(0)

        # update the plot data
        # protein
    kr_subtype_plot_protein.x_range.factors = reverseTuple(x_axis_subtype_protein)
    kr_source_dict_subtype_protein.data = dict(x=reverseTuple(x_axis_subtype_protein), y=y_axis_subtype_protein,
                                               upper=upper_bar_protein, lower=lower_bar_protein)
    # mRNA
    kr_subtype_plot_mrna.x_range.factors = reverseTuple(x_axis_subtype_mrna)
    kr_source_dict_subtype_mrna.data = dict(x=reverseTuple(x_axis_subtype_mrna), y=y_axis_subtype_mrna,
                                            upper=upper_bar_mrna, lower=lower_bar_mrna)


def me_subtype_update(attrname, old, new):
    """call back to update subtype plot for Mertins"""
    gene = new

    x_axis_subtype_protein = []
    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    x_axis_subtype_mrna = []
    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    GENE_FOR_SUBTYPES_PLOT[TICKER_INDEX] = gene  # replace the GENE_FOR_SUBTYPES_PLOT list with the new gene
    for gene_to_plot in GENE_FOR_SUBTYPES_PLOT:
        for subtype in me_subtypes:
            # protein
            x_axis_subtype_protein.append((subtype, gene_to_plot))  # get the gene name for the x axis
            try:
                y_axis_subtype_protein.append((me_protein_mean_by_subtype.loc[subtype, gene_to_plot]))
                upper_bar_protein.append(me_protein_mean_by_subtype.loc[subtype, gene_to_plot] + me_protein_sem_by_subtype.loc[subtype, gene_to_plot])
                lower_bar_protein.append(me_protein_mean_by_subtype.loc[subtype, gene_to_plot] - me_protein_sem_by_subtype.loc[subtype, gene_to_plot])
            except KeyError:
                y_axis_subtype_protein.append(0)
                upper_bar_protein.append(0)
                lower_bar_protein.append(0)

            # mRNA
            x_axis_subtype_mrna.append((subtype, gene_to_plot))
            try:
                y_axis_subtype_mrna.append((me_mrna_mean_by_subtype.loc[subtype, gene_to_plot]))
                upper_bar_mrna.append(me_mrna_mean_by_subtype.loc[subtype, gene_to_plot] + me_mrna_sem_by_subtype.loc[subtype, gene_to_plot])
                lower_bar_mrna.append(me_mrna_mean_by_subtype.loc[subtype, gene_to_plot] - me_mrna_sem_by_subtype.loc[subtype, gene_to_plot])
            except KeyError:
                y_axis_subtype_mrna.append(0)
                upper_bar_mrna.append(0)
                lower_bar_mrna.append(0)

    # update the plot data
    # protein
    me_subtype_plot_protein.x_range.factors = reverseTuple(x_axis_subtype_protein)
    me_source_dict_subtype_protein.data = dict(x=reverseTuple(x_axis_subtype_protein), y=y_axis_subtype_protein,
                                               upper=upper_bar_protein, lower=lower_bar_protein)
    # mRNA
    me_subtype_plot_mrna.x_range.factors = reverseTuple(x_axis_subtype_mrna)
    me_source_dict_subtype_mrna.data = dict(x=reverseTuple(x_axis_subtype_mrna), y=y_axis_subtype_mrna,
                                            upper=upper_bar_mrna, lower=lower_bar_mrna)


def update(attrname, old, new):
    """summary of call back functions"""
    correlation_update(attrname, old, new)
    jo_subtype_update(attrname, old, new)
    kr_subtype_update(attrname, old, new)
    me_subtype_update(attrname, old, new)


def update_mrna_correlation(attrname, old, new):
    """ function that responds to correlation table text box change"""
    Gene = new
    print("Step 1: get new gene name")

    # Johansson
    try:
        new_jo_mrna_without_input = jo_mrna_for_cor.drop(Gene) # drop the newly input gene and drop out NA
        new_jo_mrna_without_input.dropna(inplace=True)

        new_jo_mrna_cor_df = []
        for i in range(len(new_jo_mrna_without_input)):
            new_jo_mrna_cor_df.append([new_jo_mrna_without_input.index[i],
                                       round(sc.pearsonr(jo_mrna_for_cor.loc[Gene], new_jo_mrna_without_input.iloc[i])[0], 4),
                                       sc.pearsonr(jo_mrna_for_cor.loc[Gene], new_jo_mrna_without_input.iloc[i])[1]])

        new_jo_mrna_cor_df = pd.DataFrame(new_jo_mrna_cor_df, columns=['Gene', 'r', 'p']) # turn the list into df and give column labels
        new_jo_mrna_cor_df = new_jo_mrna_cor_df.sort_values('p', ascending=True) # sort the df by p values, starting from the smallest
        new_jo_top_gene_cor = new_jo_mrna_cor_df.head(100) # only display top 100 genes
    except KeyError:
        new_jo_top_gene_cor = pd.DataFrame({"Gene" : [f"No data for {Gene}"],
                                            "r": ["NA"],
                                            "p": ["NA"]})

    # Krug
    try:
        new_kr_mrna_without_input = kr_mrna_for_cor.drop(Gene) # drop the newly input gene and drop out NA
        new_kr_mrna_without_input.dropna(inplace=True)

        new_kr_mrna_cor_df = []
        for i in range(len(new_kr_mrna_without_input)):
            new_kr_mrna_cor_df.append([new_kr_mrna_without_input.index[i],
                                       round(sc.pearsonr(kr_mrna_for_cor.loc[Gene], new_kr_mrna_without_input.iloc[i])[0], 4),
                                       sc.pearsonr(kr_mrna_for_cor.loc[Gene], new_kr_mrna_without_input.iloc[i])[1]])

        new_kr_mrna_cor_df = pd.DataFrame(new_kr_mrna_cor_df, columns=['Gene', 'r', 'p']) # turn the list into df and give column labels
        new_kr_mrna_cor_df = new_kr_mrna_cor_df.sort_values('p', ascending=True) # sort the df by p values, starting from the smallest
        new_kr_top_gene_cor = new_kr_mrna_cor_df.head(100) # only display top 100 genes
    except KeyError:
        new_kr_top_gene_cor = pd.DataFrame({"Gene" : [f"No data for {Gene}"],
                                            "r": ["NA"],
                                            "p": ["NA"]})

    # Mertins
    try:
        new_me_mrna_without_input = me_mrna_for_cor.drop(Gene) # drop the newly input gene and drop out NA
        new_me_mrna_without_input.dropna(inplace=True)

        new_me_mrna_cor_df = []
        for i in range(len(new_me_mrna_without_input)):
            new_me_mrna_cor_df.append([new_me_mrna_without_input.index[i],
                                       round(sc.pearsonr(me_mrna_for_cor.loc[Gene], new_me_mrna_without_input.iloc[i])[0], 4),
                                       sc.pearsonr(me_mrna_for_cor.loc[Gene], new_me_mrna_without_input.iloc[i])[1]])

        new_me_mrna_cor_df = pd.DataFrame(new_me_mrna_cor_df, columns=['Gene', 'r', 'p']) # turn the list into df and give column labels
        new_me_mrna_cor_df = new_me_mrna_cor_df.sort_values('p', ascending=True) # sort the df by p values, starting from the smallest
        new_me_top_gene_cor = new_me_mrna_cor_df.head(100) # only display top 100 genes
    except:
        new_me_top_gene_cor = pd.DataFrame({"Gene" : [f"No data for {Gene}"],
                                            "r": ["NA"],
                                            "p": ["NA"]})

    print("Step 2: put together top 100 genes")

    jo_mrna_gene_cor_source.data = new_jo_top_gene_cor # update column data source
    jo_mrna_cor_data_table.source = jo_mrna_gene_cor_source # update table

    kr_mrna_gene_cor_source.data = new_kr_top_gene_cor  # update column data source
    kr_mrna_cor_data_table.source = kr_mrna_gene_cor_source  # update table

    me_mrna_gene_cor_source.data = new_me_top_gene_cor  # update column data source
    me_mrna_cor_data_table.source = me_mrna_gene_cor_source  # update table

    mRNA_table_gene_col_name.title = f'Genes correlated with {Gene}' # update plot title
    print("Step 3: change table")

def update_protein_correlation(attrname, old, new):
    """ function that responds to correlation table  text box change"""
    Gene = new
    print("Step 1: get new gene name")

    # Johansson
    try:
        new_jo_pro_without_input = jo_protein_for_cor.drop(Gene) # drop the newly input gene and drop out NA
        new_jo_pro_without_input.dropna(inplace=True)

        new_jo_pro_cor_df = []
        for i in range(len(new_jo_pro_without_input)):
            new_jo_pro_cor_df.append([new_jo_pro_without_input.index[i],
                                      round(sc.pearsonr(jo_protein_for_cor.loc[Gene], new_jo_pro_without_input.iloc[i])[0], 4),
                                      sc.pearsonr(jo_protein_for_cor.loc[Gene], new_jo_pro_without_input.iloc[i])[1]])

        new_jo_pro_cor_df = pd.DataFrame(new_jo_pro_cor_df, columns=['Gene', 'r', 'p']) # turn the list into df and give column labels
        new_jo_pro_cor_df = new_jo_pro_cor_df.sort_values('p', ascending=True) # sort the df by p values, starting from the smallest
        new_jo_pro_top_gene_cor = new_jo_pro_cor_df.head(100) # only display top 100 genes
    except KeyError:
        new_jo_pro_top_gene_cor = pd.DataFrame({"Gene" : [f"No data for {Gene}"],
                                            "r": ["NA"],
                                            "p": ["NA"]})

    # Krug
    try:
        new_kr_pro_without_input = kr_protein_for_cor.drop(Gene) # drop the newly input gene and drop out NA
        new_kr_pro_without_input.dropna(inplace=True)

        new_kr_pro_cor_df = []
        for i in range(len(new_kr_pro_without_input)):
            new_kr_pro_cor_df.append([new_kr_pro_without_input.index[i],
                                      round(sc.pearsonr(kr_protein_for_cor.loc[Gene], new_kr_pro_without_input.iloc[i])[0], 4),
                                      sc.pearsonr(kr_protein_for_cor.loc[Gene], new_kr_pro_without_input.iloc[i])[1]])

        new_kr_pro_cor_df = pd.DataFrame(new_kr_pro_cor_df, columns=['Gene', 'r', 'p']) # turn the list into df and give column labels
        new_kr_pro_cor_df = new_kr_pro_cor_df.sort_values('p', ascending=True) # sort the df by p values, starting from the smallest
        new_kr_pro_top_gene_cor = new_kr_pro_cor_df.head(100) # only display top 100 genes
    except KeyError:
        new_kr_pro_top_gene_cor = pd.DataFrame({"Gene" : [f"No data for {Gene}"],
                                            "r": ["NA"],
                                            "p": ["NA"]})

    # Mertins
    try:
        new_me_pro_without_input = me_protein_for_cor.drop(Gene) # drop the newly input gene and drop out NA
        new_me_pro_without_input.dropna(inplace=True)

        new_me_pro_cor_df = []
        for i in range(len(new_me_pro_without_input)):
            new_me_pro_cor_df.append([new_me_pro_without_input.index[i],
                                      round(sc.pearsonr(me_protein_for_cor.loc[Gene], new_me_pro_without_input.iloc[i])[0], 4),
                                      sc.pearsonr(me_protein_for_cor.loc[Gene], new_me_pro_without_input.iloc[i])[1]])

        new_me_pro_cor_df = pd.DataFrame(new_me_pro_cor_df, columns=['Gene', 'r', 'p']) # turn the list into df and give column labels
        new_me_pro_cor_df = new_me_pro_cor_df.sort_values('p', ascending=True) # sort the df by p values, starting from the smallest
        new_me_pro_top_gene_cor = new_me_pro_cor_df.head(100) # only display top 100 genes
    except:
        new_me_pro_top_gene_cor = pd.DataFrame({"Gene" : [f"No data for {Gene}"],
                                            "r": ["NA"],
                                            "p": ["NA"]})

    print("Step 2: put together top 100 genes")

    jo_pro_gene_cor_source.data = new_jo_pro_top_gene_cor # update column data source
    jo_pro_cor_data_table.source = jo_pro_gene_cor_source # update table

    kr_pro_gene_cor_source.data = new_kr_pro_top_gene_cor  # update column data source
    kr_pro_cor_data_table.source = kr_pro_gene_cor_source  # update table

    me_pro_gene_cor_source.data = new_me_pro_top_gene_cor  # update column data source
    me_pro_cor_data_table.source = me_pro_gene_cor_source  # update table

    protein_table_gene_col_name.title = f'Genes correlated with {Gene}' # update plot title
    print("Step 3: change table")


################################# Row 1: Protein Complex Subunit Correlation Plot ######################################
LINE_PLOT_WIDTH = 1000
LINE_PLOT_HEIGHT = 200

# Initializing the figure object
jo_plot_p = figure(x_range=FactorRange(*tuple(zip(jo_df_protein['Pam50'], jo_df_protein.index))), title='',
                   x_axis_label='', y_axis_label='protein z-score',
                   plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)
jo_plot_m = figure(x_range=FactorRange(*tuple(zip(jo_df_mrna['Pam50'], jo_df_mrna.index))), title='', x_axis_label='',
                   y_axis_label='mRNA z-score',
                   plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)

kr_plot_p = figure(x_range=FactorRange(*tuple(zip(kr_df_protein['PAM50'], kr_df_protein.index))), title='',
                   x_axis_label='', y_axis_label='protein z-score',
                   plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)
kr_plot_m = figure(x_range=FactorRange(*tuple(zip(kr_df_mrna['PAM50'], kr_df_mrna.index))), title='', x_axis_label='',
                   y_axis_label='mRNA z-score',
                   plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)

me_plot_p = figure(x_range=FactorRange(*tuple(zip(me_df_protein['PAM50'], me_df_protein.index))), title='',
                   x_axis_label='', y_axis_label='protein z-score',
                   plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)
me_plot_m = figure(x_range=FactorRange(*tuple(zip(me_df_mrna['PAM50'], me_df_mrna.index))), title='', x_axis_label='',
                   y_axis_label='mRNA z-score',
                   plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)
subtype = 'subtype'
protein_data = 'protein_data'
mRNA_data = 'mRNA_data'
gene = 'gene'

# Johansson figure
for j in range(4):
    jo_plot_p.line(subtype, protein_data, source=jo_LineSourceList[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    jo_plot_p.circle(subtype, protein_data, source=jo_LineSourceList[j], color=GENE_COLORS[j], size=4, legend_label=GENE_NUMBER[j])
    jo_plot_m.line(subtype, mRNA_data, source=jo_LineSourceList[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    jo_plot_m.circle(subtype, mRNA_data, source=jo_LineSourceList[j], color=GENE_COLORS[j], size=4, legend_label=GENE_NUMBER[j])
    TICKER[j] = AutocompleteInput(title=GENE_NUMBER[j], value=INITIAL_GENE[j], width=100, width_policy='auto',
                                  min_characters=3, completions=nix(INITIAL_GENE[j], TICKER_GENELIST), case_sensitive=False)
    TICKER[j].on_change('value', TICKER_FUNCTION_LIST[j])
    TICKER[j].on_change('value', update)
jo_plot_p = StylePlots(jo_plot_p, PlotID='Correlation')
jo_plot_m = StylePlots(jo_plot_m, PlotID='Correlation')

# Krug figure
for j in range(4):
    kr_plot_p.line(subtype, protein_data, source=kr_LineSourceList[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    kr_plot_p.circle(subtype, protein_data, source=kr_LineSourceList[j], color=GENE_COLORS[j], size=4, legend_label=GENE_NUMBER[j])
    kr_plot_m.line(subtype, mRNA_data, source=kr_LineSourceList[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    kr_plot_m.circle(subtype, mRNA_data, source=kr_LineSourceList[j], color=GENE_COLORS[j], size=4, legend_label=GENE_NUMBER[j])
    TICKER[j] = AutocompleteInput(title=GENE_NUMBER[j], value=INITIAL_GENE[j], width=100, width_policy='auto',
                                  min_characters=3, completions=nix(INITIAL_GENE[j], TICKER_GENELIST), case_sensitive=False)
    TICKER[j].on_change('value', TICKER_FUNCTION_LIST[j])
    TICKER[j].on_change('value', update)
kr_plot_p = StylePlots(kr_plot_p, PlotID='Correlation')
kr_plot_m = StylePlots(kr_plot_m, PlotID='Correlation')

# Mertins figure
for j in range(4):
    me_plot_p.line(subtype, protein_data, source=me_LineSourceList[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    me_plot_p.circle(subtype, protein_data, source=me_LineSourceList[j], color=GENE_COLORS[j], size=4, legend_label=GENE_NUMBER[j])
    me_plot_m.line(subtype, mRNA_data, source=me_LineSourceList[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    me_plot_m.circle(subtype, mRNA_data, source=me_LineSourceList[j], color=GENE_COLORS[j], size=4, legend_label=GENE_NUMBER[j])
    TICKER[j] = AutocompleteInput(title=GENE_NUMBER[j], value=INITIAL_GENE[j], width=100, width_policy='auto',
                                  min_characters=3, completions=nix(INITIAL_GENE[j], TICKER_GENELIST), case_sensitive=False)
    TICKER[j].on_change('value', TICKER_FUNCTION_LIST[j])
    TICKER[j].on_change('value', update)
me_plot_p = StylePlots(me_plot_p, PlotID='Correlation')
me_plot_m = StylePlots(me_plot_m, PlotID='Correlation')

# Put protein plot and mRNA plot into the same column, with a spacer separating them
line_plot_jo_layout = layout(column(jo_plot_p, Spacer(height=30), jo_plot_m))
line_plot_kr_layout = layout(column(kr_plot_p, Spacer(height=30), kr_plot_m))
line_plot_me_layout = layout(column(me_plot_p, Spacer(height=30), me_plot_m))

# Stack the 3 layouts into tabs
line_plot_tab = Tabs(tabs=[Panel(child=line_plot_jo_layout, title="Johansson"),
                     Panel(child=line_plot_kr_layout, title="Krug"),
                     Panel(child=line_plot_me_layout, title="Mertins")])

#Download data
button1 = Button(label="Download Gene 1 Data", button_type="success", width = 150)
button1.js_on_event("button_click",CustomJS(args=dict(source=jo_LineSourceList[0]),code=open(join(dirname(__file__),"download_javascript/jo_download.js")).read()))
button1.js_on_event("button_click",CustomJS(args=dict(source=kr_LineSourceList[0]),code=open(join(dirname(__file__),"download_javascript/kr_download.js")).read()))
button1.js_on_event("button_click",CustomJS(args=dict(source=me_LineSourceList[0]),code=open(join(dirname(__file__),"download_javascript/me_download.js")).read()))
####
button2 = Button(label="Download Gene 2 Data", button_type="success", width = 150)
button2.js_on_event("button_click", CustomJS(args=dict(source=jo_LineSourceList[1]),code=open(join(dirname(__file__),"download_javascript/jo_download.js")).read()))
button2.js_on_event("button_click", CustomJS(args=dict(source=kr_LineSourceList[1]),code=open(join(dirname(__file__),"download_javascript/kr_download.js")).read()))
button2.js_on_event("button_click", CustomJS(args=dict(source=me_LineSourceList[1]),code=open(join(dirname(__file__),"download_javascript/me_download.js")).read()))
####
button3 = Button(label="Download Gene 3 Data", button_type="success", width = 150)
button3.js_on_event("button_click", CustomJS(args=dict(source=jo_LineSourceList[2]),code=open(join(dirname(__file__),"download_javascript/jo_download.js")).read()))
button3.js_on_event("button_click", CustomJS(args=dict(source=kr_LineSourceList[2]),code=open(join(dirname(__file__),"download_javascript/kr_download.js")).read()))
button3.js_on_event("button_click", CustomJS(args=dict(source=me_LineSourceList[2]),code=open(join(dirname(__file__),"download_javascript/me_download.js")).read()))
####
button4 = Button(label="Download Gene 4 Data", button_type="success", width = 150)
button4.js_on_event("button_click", CustomJS(args=dict(source=jo_LineSourceList[3]),code=open(join(dirname(__file__),"download_javascript/jo_download.js")).read()))
button4.js_on_event("button_click", CustomJS(args=dict(source=kr_LineSourceList[3]),code=open(join(dirname(__file__), "download_javascript/kr_download.js")).read()))
button4.js_on_event("button_click", CustomJS(args=dict(source=me_LineSourceList[3]),code=open(join(dirname(__file__),"download_javascript/me_download.js")).read()))

################################# Row 2: mRNA-Protein Correlation Scatter Plot #########################################
mRNA_PRO_PLOT_WIDTH = 515
mRNA_PRO_PLOT_HEIGHT = 275

# Johansson figure
jo_plot_mRNA_prot1 = figure(title=GENE_NUMBER[0], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
jo_plot_mRNA_prot1.scatter('x1', 'y1', source=jo_cor2_source,
                           color=factor_cmap('type1', jo_subtype_colors, jo_subtypes), legend_group='type1')

jo_plot_mRNA_prot2 = figure(title=GENE_NUMBER[1], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
jo_plot_mRNA_prot2.scatter('x2', 'y2', source=jo_cor2_source,
                           color=factor_cmap('type2', jo_subtype_colors, jo_subtypes), legend_group='type2')

jo_plot_mRNA_prot3 = figure(title=GENE_NUMBER[2], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
jo_plot_mRNA_prot3.scatter('x3', 'y3', source=jo_cor2_source,
                           color=factor_cmap('type3', jo_subtype_colors, jo_subtypes), legend_group='type3')

jo_plot_mRNA_prot4 = figure(title=GENE_NUMBER[3], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
jo_plot_mRNA_prot4.scatter('x4', 'y4', source=jo_cor2_source,
                           color=factor_cmap('type4', jo_subtype_colors, jo_subtypes), legend_group='type4')

jo_plot_mRNA_prot = [jo_plot_mRNA_prot1, jo_plot_mRNA_prot2, jo_plot_mRNA_prot3, jo_plot_mRNA_prot4]

for i in jo_plot_mRNA_prot:
    StylePlots(i, PlotID='mRNA-Prot')

# Krug figure
kr_plot_mRNA_prot1 = figure(title=GENE_NUMBER[0], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
kr_plot_mRNA_prot1.scatter('x1', 'y1', source=kr_cor2_source,
                           color=factor_cmap('type1', kr_subtype_colors, kr_subtypes), legend_group='type1')

kr_plot_mRNA_prot2 = figure(title=GENE_NUMBER[1], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
kr_plot_mRNA_prot2.scatter('x2', 'y2', source=kr_cor2_source,
                           color=factor_cmap('type2', kr_subtype_colors, kr_subtypes), legend_group='type2')

kr_plot_mRNA_prot3 = figure(title=GENE_NUMBER[2], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
kr_plot_mRNA_prot3.scatter('x3', 'y3', source=kr_cor2_source,
                           color=factor_cmap('type3', kr_subtype_colors, kr_subtypes), legend_group='type3')

kr_plot_mRNA_prot4 = figure(title=GENE_NUMBER[3], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
kr_plot_mRNA_prot4.scatter('x4', 'y4', source=kr_cor2_source,
                           color=factor_cmap('type4', kr_subtype_colors, kr_subtypes), legend_group='type4')

kr_plot_mRNA_prot = [kr_plot_mRNA_prot1, kr_plot_mRNA_prot2, kr_plot_mRNA_prot3, kr_plot_mRNA_prot4]

for i in kr_plot_mRNA_prot:
    StylePlots(i, PlotID='mRNA-Prot')

# Mertins figure
me_plot_mRNA_prot1 = figure(title=GENE_NUMBER[0], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
me_plot_mRNA_prot1.scatter('x1', 'y1', source=me_cor2_source,
                           color=factor_cmap('type1', me_subtype_colors, me_subtypes), legend_group='type1')

me_plot_mRNA_prot2 = figure(title=GENE_NUMBER[1], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
me_plot_mRNA_prot2.scatter('x2', 'y2', source=me_cor2_source,
                           color=factor_cmap('type2', me_subtype_colors, me_subtypes), legend_group='type2')

me_plot_mRNA_prot3 = figure(title=GENE_NUMBER[2], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
me_plot_mRNA_prot3.scatter('x3', 'y3', source=me_cor2_source,
                           color=factor_cmap('type3', me_subtype_colors, me_subtypes), legend_group='type3')

me_plot_mRNA_prot4 = figure(title=GENE_NUMBER[3], x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                            plot_width=mRNA_PRO_PLOT_WIDTH, plot_height=mRNA_PRO_PLOT_HEIGHT)
me_plot_mRNA_prot4.scatter('x4', 'y4', source=me_cor2_source,
                           color=factor_cmap('type4', me_subtype_colors, me_subtypes), legend_group='type4')

me_plot_mRNA_prot = [me_plot_mRNA_prot1, me_plot_mRNA_prot2, me_plot_mRNA_prot3, me_plot_mRNA_prot4]

for i in me_plot_mRNA_prot:
    StylePlots(i, PlotID='mRNA-Prot')

# put gene1 and gene2 mRNA-protein correlation plots in the same row, with a spacer between
# put gene3 and gene4 mRNA-protein correlation plots in the row below, with a spacer between.
# spacer to divide the rows
scatter_plot_jo_layout = layout(column(row(jo_plot_mRNA_prot[0],jo_plot_mRNA_prot[1]),
                                       Spacer(height=30),
                                       row(jo_plot_mRNA_prot[2],jo_plot_mRNA_prot[3])))
scatter_plot_kr_layout = layout(column(row(kr_plot_mRNA_prot[0],kr_plot_mRNA_prot[1]),
                                       Spacer(height=30),
                                       row(kr_plot_mRNA_prot[2],kr_plot_mRNA_prot[3])))
scatter_plot_me_layout = layout(column(row(me_plot_mRNA_prot[0],me_plot_mRNA_prot[1]),
                                       Spacer(height=30),
                                       row(me_plot_mRNA_prot[2],me_plot_mRNA_prot[3])))
# stack 3 layouts into tab
scatter_plot_tab = Tabs(tabs=[Panel(child=scatter_plot_jo_layout, title="Johansson"),
                              Panel(child=scatter_plot_kr_layout, title="Krug"),
                              Panel(child=scatter_plot_me_layout, title="Mertins")])

########################################  Row 3: Subtype Mean and SEM plot #############################################
# set aesthetic values
jo_plot_subtype_p = StylePlots(jo_subtype_plot_protein, PlotID='Subtype')
jo_plot_subtype_m = StylePlots(jo_subtype_plot_mrna, PlotID='Subtype')

kr_plot_subtype_p = StylePlots(kr_subtype_plot_protein, PlotID='Subtype')
kr_plot_subtype_m = StylePlots(kr_subtype_plot_mrna, PlotID='Subtype')

me_plot_subtype_p = StylePlots(me_subtype_plot_protein, PlotID='Subtype')
me_plot_subtype_m = StylePlots(me_subtype_plot_mrna, PlotID='Subtype')

# Put the subtype plots for protein and mRNA into the same row, each plot takes a column, with a spacer spearating them
subtype_plot_jo_layout = layout(row(column(jo_plot_subtype_p), column(jo_plot_subtype_m)))
subtype_plot_kr_layout = layout(row(column(kr_plot_subtype_p), column(kr_plot_subtype_m)))
subtype_plot_me_layout = layout(row(column(me_plot_subtype_p), column(me_plot_subtype_m)))

# Stack the 3 layouts into tabs
subtype_plot_tab = Tabs(tabs=[Panel(child=subtype_plot_jo_layout, title="Johansson"),
                              Panel(child=subtype_plot_kr_layout, title="Krug"),
                              Panel(child=subtype_plot_me_layout, title="Mertins")])

################################################### Row 4: GxG Correlation #############################################
# set up unique genes gathered from 3 datasets
mrna_cor_gene_list = list(jo_mrna_for_cor.index) + list(kr_mrna_for_cor.index) + list(me_mrna_for_cor.index)
mrna_cor_gene_list = np.ndarray.tolist(np.unique(mrna_cor_gene_list))

protein_cor_gene_list = list(jo_protein_for_cor.index) + list(kr_protein_for_cor.index) + list(me_protein_for_cor.index)
protein_cor_gene_list = np.ndarray.tolist(np.unique(protein_cor_gene_list))

######## Johansson - mRNA
# drop the gene in the textbox, so the input gene can be correlated with the rest of the genes.
jo_mrna_without_input = jo_mrna_for_cor.drop(Gene)
jo_mrna_without_input = jo_mrna_without_input.dropna()

# create a list to collect gene name, r value and p value that the input gene is correlated with
jo_mrna_cor = []
for i in range(len(jo_mrna_without_input)):
    jo_mrna_cor.append([jo_mrna_without_input.index[i],
                        round(sc.pearsonr(jo_mrna_for_cor.loc[Gene], jo_mrna_without_input.iloc[i])[0], 4),
                        sc.pearsonr(jo_mrna_for_cor.loc[Gene], jo_mrna_without_input.iloc[i])[1]])

jo_mrna_cor = pd.DataFrame(jo_mrna_cor, columns = ['Gene', 'r', 'p'])
jo_mrna_cor = jo_mrna_cor.sort_values('p', ascending=True) # sort the gene by smallest p value
jo_mrna_top_gene_cor = jo_mrna_cor.head(100) # retrieve the top 100 most significant correlated genes

######## Krug - mRNA
# drop the gene in the textbox, so the input gene can be correlated with the rest of the genes.
kr_mrna_without_input = kr_mrna_for_cor.drop(Gene)
kr_mrna_without_input = kr_mrna_without_input.dropna()

# create a list to collect gene name, r value and p value that the input gene is correlated with
kr_mrna_cor = []
for i in range(len(kr_mrna_without_input)):
    kr_mrna_cor.append([kr_mrna_without_input.index[i],
                        round(sc.pearsonr(kr_mrna_for_cor.loc[Gene], kr_mrna_without_input.iloc[i])[0], 4),
                        sc.pearsonr(kr_mrna_for_cor.loc[Gene], kr_mrna_without_input.iloc[i])[1]])

kr_mrna_cor = pd.DataFrame(kr_mrna_cor, columns = ['Gene', 'r', 'p'])
kr_mrna_cor = kr_mrna_cor.sort_values('p', ascending=True) # sort the gene by smallest p value
kr_mrna_top_gene_cor = kr_mrna_cor.head(100) # retrieve the top 100 most significant correlated genes

######## Mertins - mRNA
# drop the gene in the textbox, so the input gene can be correlated with the rest of the genes.
me_mrna_without_input = me_mrna_for_cor.drop(Gene)
me_mrna_without_input = me_mrna_without_input.dropna()

# create a list to collect gene name, r value and p value that the input gene is correlated with
me_mrna_cor = []
for i in range(len(me_mrna_without_input)):
    me_mrna_cor.append([me_mrna_without_input.index[i],
                        round(sc.pearsonr(me_mrna_for_cor.loc[Gene], me_mrna_without_input.iloc[i])[0], 4),
                        sc.pearsonr(me_mrna_for_cor.loc[Gene], me_mrna_without_input.iloc[i])[1]])

me_mrna_cor = pd.DataFrame(me_mrna_cor, columns = ['Gene', 'r', 'p'])
me_mrna_cor = me_mrna_cor.sort_values('p', ascending=True) # sort the gene by smallest p value
me_mrna_top_gene_cor = me_mrna_cor.head(100) # retrieve the top 100 most significant correlated genes

######## Johansson - protein
# drop the gene in the textbox, so the input gene can be correlated with the rest of the genes.
jo_pro_without_input = jo_protein_for_cor.drop(Gene)
jo_pro_without_input = jo_pro_without_input.dropna()

# create a list to collect gene name, r value and p value that the input gene is correlated with
jo_pro_cor = []
for i in range(len(jo_pro_without_input)):
    jo_pro_cor.append([jo_pro_without_input.index[i],
                       round(sc.pearsonr(jo_protein_for_cor.loc[Gene], jo_pro_without_input.iloc[i])[0], 4),
                       sc.pearsonr(jo_protein_for_cor.loc[Gene], jo_pro_without_input.iloc[i])[1]])

jo_pro_cor = pd.DataFrame(jo_pro_cor, columns = ['Gene', 'r', 'p'])
jo_pro_cor = jo_pro_cor.sort_values('p', ascending=True) # sort the gene by smallest p value
jo_pro_top_gene_cor = jo_pro_cor.head(100) # retrieve the top 100 most significant correlated genes

######## Krug - protein
# drop the gene in the textbox, so the input gene can be correlated with the rest of the genes.
kr_pro_without_input = kr_protein_for_cor.drop(Gene)
kr_pro_without_input = kr_pro_without_input.dropna()

# create a list to collect gene name, r value and p value that the input gene is correlated with
kr_pro_cor = []
for i in range(len(kr_pro_without_input)):
    kr_pro_cor.append([kr_pro_without_input.index[i],
                       round(sc.pearsonr(kr_protein_for_cor.loc[Gene], kr_pro_without_input.iloc[i])[0], 4),
                       sc.pearsonr(kr_protein_for_cor.loc[Gene], kr_pro_without_input.iloc[i])[1]])

kr_pro_cor = pd.DataFrame(kr_pro_cor, columns = ['Gene', 'r', 'p'])
kr_pro_cor = kr_pro_cor.sort_values('p', ascending=True) # sort the gene by smallest p value
kr_pro_top_gene_cor = kr_pro_cor.head(100) # retrieve the top 100 most significant correlated genes

######## Mertins - protein
# drop the gene in the textbox, so the input gene can be correlated with the rest of the genes.
me_pro_without_input = me_protein_for_cor.drop(Gene)
me_pro_without_input = me_pro_without_input.dropna()

# create a list to collect gene name, r value and p value that the input gene is correlated with
me_pro_cor = []
for i in range(len(me_pro_without_input)):
    me_pro_cor.append([me_pro_without_input.index[i],
                       round(sc.pearsonr(me_protein_for_cor.loc[Gene], me_pro_without_input.iloc[i])[0], 4),
                       sc.pearsonr(me_protein_for_cor.loc[Gene], me_pro_without_input.iloc[i])[1]])

me_pro_cor = pd.DataFrame(me_pro_cor, columns = ['Gene', 'r', 'p'])
me_pro_cor = me_pro_cor.sort_values('p', ascending=True) # sort the gene by smallest p value
me_pro_top_gene_cor = me_pro_cor.head(100) # retrieve the top 100 most significant correlated genes

######## Widget and aesthetic

# update when textbox input changes
mrna_correlation_textbox = AutocompleteInput(title='Gene Name for mRNA correlation table:', value=str(Gene), width=150, width_policy ='auto',
                                             min_characters=3, completions = mrna_cor_gene_list, case_sensitive=False)
mrna_correlation_textbox.on_change('value', update_mrna_correlation)

protein_correlation_textbox = AutocompleteInput(title='Gene Name for protein correlation table:', value=str(Gene), width=150, width_policy ='auto',
                                             min_characters=3, completions = protein_cor_gene_list, case_sensitive=False)
protein_correlation_textbox.on_change('value', update_protein_correlation)

# aesthetic
jo_mrna_gene_cor_source = ColumnDataSource(data=jo_mrna_top_gene_cor) # set up columndatasource
kr_mrna_gene_cor_source = ColumnDataSource(data=kr_mrna_top_gene_cor)
me_mrna_gene_cor_source = ColumnDataSource(data=me_mrna_top_gene_cor)

jo_pro_gene_cor_source = ColumnDataSource(data=jo_pro_top_gene_cor) # set up columndatasource
kr_pro_gene_cor_source = ColumnDataSource(data=kr_pro_top_gene_cor)
me_pro_gene_cor_source = ColumnDataSource(data=me_pro_top_gene_cor)

mRNA_table_gene_col_name = TableColumn(field='Gene', title=f'Genes correlated with {mrna_correlation_textbox.value}')
protein_table_gene_col_name = TableColumn(field='Gene', title=f'Genes correlated with {protein_correlation_textbox.value}')

mRNA_table_columns = [mRNA_table_gene_col_name,  # set up mRNA_table_columns
                      TableColumn(field='r', title='Coefficient'),
                      TableColumn(field='p', title='p value')]

protein_table_columns = [protein_table_gene_col_name,  # set up mRNA_table_columns
                         TableColumn(field='r', title='Coefficient'),
                         TableColumn(field='p', title='p value')]

# create table widget by assembling columndatasource and mRNA_table_columns
TABLE_WIDTH = 500
TABLE_HEIGHT = 300

jo_mrna_cor_data_table = DataTable(source=jo_mrna_gene_cor_source, columns = mRNA_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)
kr_mrna_cor_data_table = DataTable(source=kr_mrna_gene_cor_source, columns = mRNA_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)
me_mrna_cor_data_table = DataTable(source=me_mrna_gene_cor_source, columns = mRNA_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)

jo_pro_cor_data_table = DataTable(source=jo_pro_gene_cor_source, columns = protein_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)
kr_pro_cor_data_table = DataTable(source=kr_pro_gene_cor_source, columns = protein_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)
me_pro_cor_data_table = DataTable(source=me_pro_gene_cor_source, columns = protein_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)

# combine tables into tabs
cor_pro_table_tab = Tabs(tabs=[Panel(child=jo_pro_cor_data_table, title ="Johansson"),
                           Panel(child=kr_pro_cor_data_table, title = "Krug"),
                           Panel(child=me_pro_cor_data_table, title = "Mertins")])

cor_mrna_table_tab = Tabs(tabs=[Panel(child=jo_mrna_cor_data_table, title ="Johansson"),
                           Panel(child=kr_mrna_cor_data_table, title = "Krug"),
                           Panel(child=me_mrna_cor_data_table, title = "Mertins")])

# button and call back
mrna_cor_data_download_button = Button(label="Download mRNA Correlation Table", button_type="success", width = 150)
mrna_cor_data_download_button.js_on_event("button_click", CustomJS(args=dict(source=jo_mrna_gene_cor_source), code=open(join(dirname(__file__),
                                                                                                                             "download_javascript/jo_cor_table_download.js")).read()))
mrna_cor_data_download_button.js_on_event("button_click", CustomJS(args=dict(source=kr_mrna_gene_cor_source), code=open(join(dirname(__file__),
                                                                                                                             "download_javascript/kr_cor_table_download.js")).read()))
mrna_cor_data_download_button.js_on_event("button_click", CustomJS(args=dict(source=me_mrna_gene_cor_source), code=open(join(dirname(__file__),
                                                                                                                             "download_javascript/me_cor_table_download.js")).read()))

pro_cor_data_download_button = Button(label="Download Protein Correlation Table", button_type="success", width = 150)
pro_cor_data_download_button.js_on_event("button_click", CustomJS(args=dict(source=jo_pro_gene_cor_source), code=open(join(dirname(__file__),
                                                                                                                             "download_javascript/jo_cor_table_download.js")).read()))
pro_cor_data_download_button.js_on_event("button_click", CustomJS(args=dict(source=kr_pro_gene_cor_source), code=open(join(dirname(__file__),
                                                                                                                             "download_javascript/kr_cor_table_download.js")).read()))
pro_cor_data_download_button.js_on_event("button_click", CustomJS(args=dict(source=me_pro_gene_cor_source), code=open(join(dirname(__file__),
                                                                                                                             "download_javascript/me_cor_table_download.js")).read()))


###################################################### Run the application #############################################
ToolTitleText = "Breast Cancer Data Portal"
ToolTitleDiv = Div(text=ToolTitleText,
                   style={'font-size': '200%', 'color': 'black', 'font-style': 'normal', 'font-weight': 'bold'},
                   width=1000)


KrugInfo = '''The BRCA cohort from Krug et al. contains 122 treatment-naive primary breast cancer tumors. Sample subtypes include 29 basal-like, 14 Her2-enriched, 57 LumA, 17 LumB and 5 normal-like samples.
             The dataset provided identifications of fully quantified 10734 gene transcripts and 7583 proteins.'''
JohanssonInfo = '''The Oslo cohort from Johansson et al. has 45 breast cancer tumors including 5 subtypes (basal-like, Her2-enriched, LumA, LumB and normal-like). Each subtype contains 9 samples.
             The dataset provided identifications of fully quantified 22581 gene transcripts and 9995 proteins.'''
MertinsInfo = '''The TCGA-BRCA cohort from Mertins et al. contains 77 breast cancer tumors. Subtypes include 18 basal-like, 12 Her2-enriched, 23 LumA and 24 LumB samples.
             The dataset provided identifications of fully quantified 20530 gene transcripts, 6359 proteins.'''
Background_InfoDiv = column(Div(text="Background Information:",
                                style={'font-size': '120%', 'color': 'black', 'font-style': 'italic',
                                       'font-weight': 'bold'}),
                            Div(text=JohanssonInfo,
                                style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}),
                            Div(text=KrugInfo,
                                style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}),
                            Div(text=MertinsInfo,
                                style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}),
                            width=1000)


DescriptiveTextLine = "Abundances of proteins that are part of the same complex are tightly correlated across breast tumors. " \
                             "However, this does not appear to be the case for the corresponding mRNA transcripts. " \
                             "The plots are initialized showing abundances of structural proteins of complex I of the electron transport chain. " \
                             "If you wish to hide the curve for a specific gene, click on the corresponding legend entries."
CorrelationTextDiv = column(Div(text="Protein Complex Subunit Correlation:",
                                style={'font-size': '120%', 'color': 'black', 'font-style': 'italic',
                                       'font-weight': 'bold'}),
                            Div(text=DescriptiveTextLine,
                                style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}),
                            width=1000)


DescriptiveTextScatter = "Protein and transcript abundances often are not correlated; indicating substantial utilization of post-transcriptional regulatory mechanisms."
mRNA_Prot_TextDiv = column(Div(text="mRNA-Protein Correlation:",
                               style={'font-size': '120%', 'color': 'black', 'font-style': 'italic',
                                      'font-weight': 'bold'}),
                           Div(text=DescriptiveTextScatter,
                               style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}),
                           width=1000)


DescriptiveTextSubtypes = "Breast cancer subtypes are defined by their gene expression profiles. " \
                          "The plots are initialized above showing abundances of ER (ESR1), PR (PGR), HER2 (ERBB2), and KI-67 (MKI67); " \
                          "four immunohistochemical markers commonly used in the clinic. " \
                          "Values are means +/- standard error of the mean (SEM). " \
                          "If data are not available for a gene, values will appear as all zeros with no SEMs."
SubtypesTextDiv = column(Div(text="Protein and mRNA Expression by Breast Cancer PAM50 Subtype:",
                             style={'font-size': '120%', 'color': 'black', 'font-style': 'italic',
                                    'font-weight': 'bold'}),
                         Div(text=DescriptiveTextSubtypes,
                             style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}),
                         width=1000)


Protein_Table_TitleText = "Correlation Table - Protein"
Mrna_Table_TitleText = "Correlation Table - mRNA"
Protein_Table_TitleDiv = column(Div(text=Protein_Table_TitleText,
                                    style={'font-size': '120%', 'color': 'black', 'font-weight': 'bold'}))
mrna_Table_TitleDiv = column(Div(text=Mrna_Table_TitleText,
                                 style={'font-size': '120%', 'color': 'black', 'font-weight': 'bold'}))
DescriptiveTextCorTable = "The gene-gene correlation is calculated from protein and mRNA expression data using the pearson correlation. " \
               "The table shows the top 100 genes that are correlated with the input gene based on the p-value. " \
               "The computation roughly takes 30 seconds, please be patient."
CorTableTextDiv = column(Div(text="Gene-Gene Correlation Table:",
                             style={'font-size': '120%', 'color': 'black', 'font-style': 'italic',
                                    'font-weight': 'bold'}),
                         Div(text=DescriptiveTextCorTable,
                             style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}),
                         width=1000)


RowSpacer = Spacer(height=30)
PageMargin = Spacer(width=30)

InstructionsText = "To view data for a different gene, type a HGNC gene symbol in the textbox."
BlankGeneInstructionText = "If you wish to leave certain textboxes blank, please input the corresponding blank entries."
WarningText = "Please do not enter same gene/blank entry across multiple textboxes."
DownloadText = 'Click "download" to save the protein and mRNA data for corresponding genes from the three studies. ' \
              'Please use Firefox or Google Chrome for better download experience.'
RedTextDiv = column(Div(text=InstructionsText, style={'font-size': '100%', 'color': 'red', 'font-style': 'italic'},
                      width=300),
                    Div(text=BlankGeneInstructionText,
                        style={'font-size': '100%', 'color': 'red', 'font-style': 'italic'},
                        width=300),
                    Div(text=WarningText, style={'font-size': '100%', 'color': 'red', 'font-style': 'italic'},
                      width=300),
                    Div(text=DownloadText, style={'font-size': '100%', 'color': 'red', 'font-style': 'italic'},
                        width=300))

formatted_button1 = row(column(Spacer(height=19),button1))
formatted_button2 = row(column(Spacer(height=19),button2))
formatted_button3 = row(column(Spacer(height=19),button3))
formatted_button4 = row(column(Spacer(height=19),button4))
TextBox_and_buttons = column(row(TICKER[0],formatted_button1), row(TICKER[1],formatted_button2), row(TICKER[2],formatted_button3), row(TICKER[3],formatted_button4))


GeneColumn = column(TextBox_and_buttons, RedTextDiv)
Line_Plots = column(row(line_plot_tab), row(CorrelationTextDiv))
mRNA_Protein_Plots = column(row(scatter_plot_tab), row(mRNA_Prot_TextDiv))
Subtype_Plots = column(row(subtype_plot_tab), row(SubtypesTextDiv))
GxGCorTables = column(row(column(Protein_Table_TitleDiv, protein_correlation_textbox, cor_pro_table_tab, pro_cor_data_download_button),
                          Spacer(width=30),
                          column(mrna_Table_TitleDiv, mrna_correlation_textbox, cor_mrna_table_tab, mrna_cor_data_download_button)),
                      row(CorTableTextDiv))

l = layout([
    [PageMargin, ToolTitleDiv],
    [RowSpacer],
    [PageMargin, Background_InfoDiv],
    [RowSpacer],
    [PageMargin, Line_Plots, GeneColumn],
    [RowSpacer],
    [PageMargin, mRNA_Protein_Plots],
    [RowSpacer],
    [PageMargin, Subtype_Plots],
    [RowSpacer],
    [PageMargin, GxGCorTables],
    [RowSpacer]])

# Create the bokeh server application
curdoc().add_root(l)

# run line in terminal
# first, cd /Users/zhuoheng/PycharmProjects/Breast Cancer Data Portal
# bokeh serve --show myapp --session-token-expiration=1000000
