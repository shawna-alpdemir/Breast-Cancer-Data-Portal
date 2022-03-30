from bokeh.plotting import figure
from bokeh.transform import factor_cmap

from Plot_Line_Scatter_ColumnDataSource import Johansson_CDS, Krug_CDS, Mertins_CDS


#################################### protein-protein/RNA-RNA scatter plot ##############################################
# this plot is make from dict 5 and dict 6 from Plot_ColumnDataSource

# constants
SCATTER_PLOT_WIDTH = 515
SCATTER_PLOT_HEIGHT = 275
GENE1 = 'NDUFS2'
GENE2 = 'NDUFS3'
GENE3 = 'NDUFS7'
GENE4 = 'blank 4'
five_subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']  # list of subtypes
five_subtype_colors = ['#E31A1C', '#FB9A99', '#1F78B4', '#A6CEE3','#33A02C']

# dictionary key name
subtype_tuple = 'subtype, tumor'
subtype = 'subtype'
protein_data = 'protein_data'
mRNA_data = 'mRNA_data'
gene = 'gene'

# function calls
[johansson_cds, johansson_subtype_tumor_tuple, jo_unique_gene_list] = Johansson_CDS()
[krug_cds, krug_subtype_tumor_tuple, kr_unique_gene_list] = Krug_CDS()
[mertins_cds, mertins_subtype_tumor_tuple, me_unique_gene_list] = Mertins_CDS()

def Pro_Pro_Scatter_Plot():
    """protein-protein for gene 1 and gene 2 scatter plot for Johansson, Krug and Mertins"""

    jo_pro_pro_scatter_plot = figure(title='', x_axis_label=f'{GENE1} protein z-score', y_axis_label=f'{GENE2} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    jo_pro_pro_scatter_plot.scatter('x_protein_data', 'y_protein_data', source=johansson_cds[4],
                               color=factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group=subtype)

    kr_pro_pro_scatter_plot = figure(title='', x_axis_label=f'{GENE1} protein z-score', y_axis_label=f'{GENE2} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    kr_pro_pro_scatter_plot.scatter('x_protein_data', 'y_protein_data', source=krug_cds[4],
                               color=factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group=subtype)

    me_pro_pro_scatter_plot = figure(title='', x_axis_label=f'{GENE1} protein z-score', y_axis_label=f'{GENE2} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    me_pro_pro_scatter_plot.scatter('x_protein_data', 'y_protein_data', source=mertins_cds[4],
                                    color=factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]), legend_group=subtype)

    return jo_pro_pro_scatter_plot, kr_pro_pro_scatter_plot, me_pro_pro_scatter_plot

def RNA_RNA_Scatter_Plot():
    """RNA-RNA for gene 1 and gene 2 scatter plot for Johansson, Krug and Mertins"""

    jo_rna_rna_scatter_plot = figure(title='', x_axis_label=f'{GENE1} mRNA z-score', y_axis_label=f'{GENE2} mRNA z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    jo_rna_rna_scatter_plot.scatter('x_mRNA_data', 'y_mRNA_data', source=johansson_cds[5],
                                color = factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group = subtype)

    kr_rna_rna_scatter_plot = figure(title='', x_axis_label=f'{GENE1} mRNA z-score', y_axis_label=f'{GENE2} mRNA z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    kr_rna_rna_scatter_plot.scatter('x_mRNA_data', 'y_mRNA_data', source=krug_cds[5],
                                color = factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group = subtype)

    me_rna_rna_scatter_plot = figure(title='', x_axis_label=f'{GENE1} mRNA z-score', y_axis_label=f'{GENE2} mRNA z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    me_rna_rna_scatter_plot.scatter('x_mRNA_data', 'y_mRNA_data', source=mertins_cds[5],
                                color = factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]), legend_group = subtype)

    return jo_rna_rna_scatter_plot, kr_rna_rna_scatter_plot, me_rna_rna_scatter_plot

def Johansson_RNA_Pro_Scatter_Plot():
    """ mRNA-Protein scatter plot for gene 1-4 for Johansson"""
    jo_plot_mRNA_prot1 = figure(title='', x_axis_label=f'{GENE1} mRNA z-score', y_axis_label=f'{GENE1} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    jo_plot_mRNA_prot1.scatter(mRNA_data, protein_data, source=johansson_cds[0],
                               color=factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group=subtype)

    jo_plot_mRNA_prot2 = figure(title='', x_axis_label=f'{GENE2} mRNA z-score', y_axis_label=f'{GENE2} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    jo_plot_mRNA_prot2.scatter(mRNA_data, protein_data, source=johansson_cds[1],
                               color=factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group=subtype)

    jo_plot_mRNA_prot3 = figure(title='', x_axis_label=f'{GENE3} mRNA z-score', y_axis_label=f'{GENE3} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    jo_plot_mRNA_prot3.scatter(mRNA_data, protein_data, source=johansson_cds[2],
                               color=factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group=subtype)

    jo_plot_mRNA_prot4 = figure(title='', x_axis_label=f'{GENE4} mRNA z-score', y_axis_label=f'{GENE4} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    jo_plot_mRNA_prot4.scatter(mRNA_data, protein_data, source=johansson_cds[3],
                               color=factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group=subtype)

    return jo_plot_mRNA_prot1, jo_plot_mRNA_prot2, jo_plot_mRNA_prot3, jo_plot_mRNA_prot4


def Krug_RNA_Pro_Scatter_Plot():
    """ mRNA-Protein scatter plot for gene 1-4 for Krug"""

    kr_plot_mRNA_prot1 = figure(title='', x_axis_label=f'{GENE1} mRNA z-score', y_axis_label=f'{GENE1} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    kr_plot_mRNA_prot1.scatter(mRNA_data, protein_data, source=krug_cds[0],
                               color=factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group=subtype)

    kr_plot_mRNA_prot2 = figure(title='', x_axis_label=f'{GENE2} mRNA z-score', y_axis_label=f'{GENE2} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    kr_plot_mRNA_prot2.scatter(mRNA_data, protein_data, source=krug_cds[1],
                               color=factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group=subtype)

    kr_plot_mRNA_prot3 = figure(title='', x_axis_label=f'{GENE3} mRNA z-score', y_axis_label=f'{GENE3} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    kr_plot_mRNA_prot3.scatter(mRNA_data, protein_data, source=krug_cds[2],
                               color=factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group=subtype)

    kr_plot_mRNA_prot4 = figure(title='', x_axis_label=f'{GENE4} mRNA z-score', y_axis_label=f'{GENE4} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    kr_plot_mRNA_prot4.scatter(mRNA_data, protein_data, source=krug_cds[3],
                               color=factor_cmap(subtype, five_subtype_colors, five_subtypes), legend_group=subtype)

    return kr_plot_mRNA_prot1, kr_plot_mRNA_prot2, kr_plot_mRNA_prot3, kr_plot_mRNA_prot4


def Mertins_RNA_Pro_Scatter_Plot():
    """ mRNA-Protein scatter plot for gene 1-4 for Mertins"""

    me_plot_mRNA_prot1 = figure(title='', x_axis_label=f'{GENE1} mRNA z-score', y_axis_label=f'{GENE1} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    me_plot_mRNA_prot1.scatter(mRNA_data, protein_data, source=mertins_cds[0],
                               color=factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]),
                               legend_group=subtype)

    me_plot_mRNA_prot2 = figure(title='', x_axis_label=f'{GENE2} mRNA z-score', y_axis_label=f'{GENE2} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    me_plot_mRNA_prot2.scatter(mRNA_data, protein_data, source=mertins_cds[1],
                               color=factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]),
                               legend_group=subtype)

    me_plot_mRNA_prot3 = figure(title='', x_axis_label=f'{GENE3} mRNA z-score', y_axis_label=f'{GENE3} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    me_plot_mRNA_prot3.scatter(mRNA_data, protein_data, source=mertins_cds[2],
                               color=factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]),
                               legend_group=subtype)

    me_plot_mRNA_prot4 = figure(title='', x_axis_label=f'{GENE4} mRNA z-score', y_axis_label=f'{GENE4} protein z-score',
                                plot_width=SCATTER_PLOT_WIDTH, plot_height=SCATTER_PLOT_HEIGHT)
    me_plot_mRNA_prot4.scatter(mRNA_data, protein_data, source=mertins_cds[3],
                               color=factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]),
                               legend_group=subtype)

    return me_plot_mRNA_prot1, me_plot_mRNA_prot2, me_plot_mRNA_prot3, me_plot_mRNA_prot4