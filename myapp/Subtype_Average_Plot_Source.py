from bokeh.models import ColumnDataSource
from functools import lru_cache # a function that helps in reducing the execution time of the function by using memoization technique
# import Import_data as im
# import Subtype_Average_DF as su

def reverseTuple(listOfTuple):
    """For reverse the list of tuple for x-axis data [(Basal, geneA)] --> [(geneA, Basal)]"""
    return list(map(lambda tup: tup[::-1], listOfTuple))


########################################################################################################################
##################################################### Johansson ########################################################

def jo_Subtype_Avg_Plot_Source(jo_protein_mean_by_subtype, jo_mrna_mean_by_subtype,
                               jo_protein_sem_by_subtype, jo_mrna_sem_by_subtype,
                               jo_subtypes, gene_for_subtypes_plot=None):
    """ Create ColumnDataSource dictionary for plotting subtype mean and SEM purposes """
    # Create empty list to store data iterated from the for loop below

    if gene_for_subtypes_plot is None:
        gene_for_subtypes_plot = ['ESR1', 'PGR', 'ERBB2', 'MKI67']
    x_axis_subtype_protein = []
    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    LegendGroup = []

    x_axis_subtype_mrna = []
    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    for gene in gene_for_subtypes_plot:
        for subtype in jo_subtypes:
            # find the MEAN and SEM data of each gene in the gene_for_subtypes_plot list
            # gene_for_subtypes carries initial four genes, and can be modified by user-input in main.py

            # legend list
            LegendGroup.append(subtype)

            # Protein lists
            x_axis_subtype_protein.append((subtype, gene))
            y_axis_subtype_protein.append(jo_protein_mean_by_subtype.loc[subtype, gene])
            upper_bar_protein.append(
                jo_protein_mean_by_subtype.loc[subtype, gene] + jo_protein_sem_by_subtype.loc[subtype, gene])
            lower_bar_protein.append(
                jo_protein_mean_by_subtype.loc[subtype, gene] - jo_protein_sem_by_subtype.loc[subtype, gene])

            # mRNA lists
            x_axis_subtype_mrna.append((subtype, gene))
            y_axis_subtype_mrna.append(jo_mrna_mean_by_subtype.loc[subtype, gene])
            upper_bar_mrna.append(
                jo_mrna_mean_by_subtype.loc[subtype, gene] + jo_mrna_sem_by_subtype.loc[subtype, gene])
            lower_bar_mrna.append(
                jo_mrna_mean_by_subtype.loc[subtype, gene] - jo_mrna_sem_by_subtype.loc[subtype, gene])

    # create the plot data sources from the dictionaries for each gene.
    jo_source_dict_subtype_protein = ColumnDataSource(
        data={'x': reverseTuple(x_axis_subtype_protein), 'y': y_axis_subtype_protein,
              'upper': upper_bar_protein, 'lower': lower_bar_protein, 'legend_group': LegendGroup})

    jo_source_dict_subtype_mrna = ColumnDataSource(
        data={'x': reverseTuple(x_axis_subtype_mrna), 'y': y_axis_subtype_mrna,
              'upper': upper_bar_mrna, 'lower': lower_bar_mrna, 'legend_group': LegendGroup})
    # The data dictionary will look like this:
    # 'x': [(geneA, basal), (geneA, her2)....]
    # 'y': [mean for geneA x basal, mean for geneA x her2...]
    # 'upper' [mean+SEM for geneA x basal, mean+SEM for geneA x her2...]
    # 'lower' [mean-SEM for geneA x basal, mean-SEM for geneA x her2...]
    # 'legend_group = [5 subtypes]

    return jo_source_dict_subtype_protein, jo_source_dict_subtype_mrna


########################################################################################################################
##################################################### Krug #############################################################

def kr_Subtype_Avg_Plot_Source(kr_protein_mean_by_subtype, kr_mrna_mean_by_subtype,
                               kr_protein_sem_by_subtype, kr_mrna_sem_by_subtype,
                               kr_subtypes, gene_for_subtypes_plot=None):
    """ Create ColumnDataSource dictionary for plotting subtype mean and SEM purposes """
    # Create empty list to store data iterated from the for loop below
    if gene_for_subtypes_plot is None:
        gene_for_subtypes_plot = ['ESR1', 'PGR', 'ERBB2', 'MKI67']
    x_axis_subtype_protein = []
    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    LegendGroup = []

    x_axis_subtype_mrna = []
    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    for gene in gene_for_subtypes_plot:
        for subtype in kr_subtypes:
            # find the MEAN and SEM data of each gene in the gene_for_subtypes_plot list
            # gene_for_subtypes carries initial four genes, and can be modified by user-input in main.py

            # legend list
            LegendGroup.append(subtype)

            # Protein lists
            x_axis_subtype_protein.append((subtype, gene))
            y_axis_subtype_protein.append(kr_protein_mean_by_subtype.loc[subtype, gene])
            upper_bar_protein.append(
                kr_protein_mean_by_subtype.loc[subtype, gene] + kr_protein_sem_by_subtype.loc[subtype, gene])
            lower_bar_protein.append(
                kr_protein_mean_by_subtype.loc[subtype, gene] - kr_protein_sem_by_subtype.loc[subtype, gene])

            # mRNA lists
            x_axis_subtype_mrna.append((subtype, gene))
            y_axis_subtype_mrna.append(kr_mrna_mean_by_subtype.loc[subtype, gene])
            upper_bar_mrna.append(
                kr_mrna_mean_by_subtype.loc[subtype, gene] + kr_mrna_sem_by_subtype.loc[subtype, gene])
            lower_bar_mrna.append(
                kr_mrna_mean_by_subtype.loc[subtype, gene] - kr_mrna_sem_by_subtype.loc[subtype, gene])

    # create the plot data sources from the dictionaries for each gene.
    kr_source_dict_subtype_protein = ColumnDataSource(
        data={'x': reverseTuple(x_axis_subtype_protein), 'y': y_axis_subtype_protein,
              'upper': upper_bar_protein, 'lower': lower_bar_protein, 'legend_group': LegendGroup})

    kr_source_dict_subtype_mrna = ColumnDataSource(
        data={'x': reverseTuple(x_axis_subtype_mrna), 'y': y_axis_subtype_mrna,
              'upper': upper_bar_mrna, 'lower': lower_bar_mrna, 'legend_group': LegendGroup})

    # The data dictionary will look like this:
    # 'x': [(geneA, basal), (geneA, her2)....]
    # 'y': [mean for geneA x basal, mean for geneA x her2...]
    # 'upper' [mean+SEM for geneA x basal, mean+SEM for geneA x her2...]
    # 'lower' [mean-SEM for geneA x basal, mean-SEM for geneA x her2...]
    # 'legend_group = [5 subtypes]

    return kr_source_dict_subtype_protein, kr_source_dict_subtype_mrna


########################################################################################################################
##################################################### Mertins ##########################################################

def me_Subtype_Avg_Plot_Source(me_protein_mean_by_subtype, me_mrna_mean_by_subtype,
                               me_protein_sem_by_subtype, me_mrna_sem_by_subtype,
                               me_subtypes, gene_for_subtypes_plot=None):
    """ Create ColumnDataSource dictionary for plotting subtype mean and SEM purposes """
    # Create empty list to store data iterated from the for loop below
    if gene_for_subtypes_plot is None:
        gene_for_subtypes_plot = ['ESR1', 'PGR', 'ERBB2', 'MKI67']
    x_axis_subtype_protein = []
    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    LegendGroup = []

    x_axis_subtype_mrna = []
    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    for gene in gene_for_subtypes_plot:
        for subtype in me_subtypes:
            # find the MEAN and SEM data of each gene in the gene_for_subtypes_plot list
            # gene_for_subtypes carries initial four genes, and can be modified by user-input in main.py

            # legend list
            LegendGroup.append(subtype)

            # Protein lists
            x_axis_subtype_protein.append((subtype, gene))
            y_axis_subtype_protein.append(me_protein_mean_by_subtype.loc[subtype, gene])
            upper_bar_protein.append(
                me_protein_mean_by_subtype.loc[subtype, gene] + me_protein_sem_by_subtype.loc[subtype, gene])
            lower_bar_protein.append(
                me_protein_mean_by_subtype.loc[subtype, gene] - me_protein_sem_by_subtype.loc[subtype, gene])

            # mRNA lists
            x_axis_subtype_mrna.append((subtype, gene))
            y_axis_subtype_mrna.append(me_mrna_mean_by_subtype.loc[subtype, gene])
            upper_bar_mrna.append(
                me_mrna_mean_by_subtype.loc[subtype, gene] + me_mrna_sem_by_subtype.loc[subtype, gene])
            lower_bar_mrna.append(
                me_mrna_mean_by_subtype.loc[subtype, gene] - me_mrna_sem_by_subtype.loc[subtype, gene])

    # create the plot data sources from the dictionaries for each gene.
    me_source_dict_subtype_protein = ColumnDataSource(
        data={'x': reverseTuple(x_axis_subtype_protein), 'y': y_axis_subtype_protein,
              'upper': upper_bar_protein, 'lower': lower_bar_protein, 'legend_group': LegendGroup})

    me_source_dict_subtype_mrna = ColumnDataSource(
        data={'x': reverseTuple(x_axis_subtype_mrna), 'y': y_axis_subtype_mrna,
              'upper': upper_bar_mrna, 'lower': lower_bar_mrna, 'legend_group': LegendGroup})
    # The data dictionary will look like this:
    # 'x': [(geneA, basal), (geneA, her2)....]
    # 'y': [mean for geneA x basal, mean for geneA x her2...]
    # 'upper' [mean+SEM for geneA x basal, mean+SEM for geneA x her2...]
    # 'lower' [mean-SEM for geneA x basal, mean-SEM for geneA x her2...]
    # 'legend_group = [5 subtypes]

    return me_source_dict_subtype_protein, me_source_dict_subtype_mrna

# # testing
# # return jo_df_protein, jo_df_mrna, jo_genes, jo_subtypes, jo_subtype_colors
# a = im.jo_ImportData()
#
# # return jo_protein_mean_by_subtype, jo_mrna_mean_by_subtype, jo_protein_sem_by_subtype, jo_mrna_sem_by_subtype
# b=su.jo_Subtype_Avg_SEM_DFs(a[0],a[1],a[3])
#
# # return jo_source_dict_subtype_protein, jo_source_dict_subtype_mrna, jo_subtypes
# test = jo_Subtype_Avg_Plot_Source(b[0],b[1],b[2],b[3],a[3], ['ESR1','PGR','ERBB2','MKI67'])
#
# #print(test[0])
#
#
# from bokeh.models import ColumnDataSource, Whisker, FactorRange
# from bokeh.plotting import figure, show
# from bokeh.transform import factor_cmap
#
# p = figure(x_range=FactorRange(*test[0].data["x"]), width=600, height=300, title="Bar graph showing protein expression of genes",
#            y_axis_label = 'protein z-score')
#
# p.circle('x', 'y', source=test[0].data,
#                          fill_color=factor_cmap('x', palette=a[4], factors=a[3], start=1, end=2),
#                          line_color=factor_cmap('x', palette=a[4], factors=a[3], start=1, end=2),
#                          legend_group='legend_group')
#
# W_p = Whisker(source=test[0], base="x", upper="upper", lower="lower",
#               line_color=factor_cmap('x', palette=a[4], factors=a[3], start=1, end=2))
# W_p.upper_head.line_color = factor_cmap('x', palette=a[4], factors=a[3], start=1, end=2)
# W_p.lower_head.line_color = factor_cmap('x', palette=a[4], factors=a[3], start=1, end=2)
# p.add_layout(W_p)
#
# p.xgrid.grid_line_color = None
#
#
# show(p)
