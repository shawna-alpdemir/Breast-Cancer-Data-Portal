# function that initiates scatter plots

from bokeh.plotting import figure

#################################### protein-protein/RNA-RNA scatter plot ##############################################
# this plot is make from dict 5 and dict 6 from Plot_ColumnDataSource

# constants
SCATTER_PLOT_WIDTH = 230
SCATTER_PLOT_HEIGHT = 230
GENE1 = 'NDUFS2'
GENE2 = 'NDUFS3'
GENE3 = 'NDUFS7'
GENE4 = 'NDUFS8'
INITIAL_GENE = ['NDUFS2', 'NDUFS3', 'NDUFS7', 'NDUFS8']
five_subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']  # list of subtypes
five_subtype_colors = ['#E31A1C', '#FB9A99', '#1F78B4', '#A6CEE3','#33A02C']

# dictionary key name
subtype_tuple = 'subtype, tumor'
subtype = 'subtype'
protein_data = 'protein_data'
mRNA_data = 'mRNA_data'
gene = 'gene'


def Pro_Pro_Scatter_Plot():
    """protein-protein for gene 1 and gene 2 scatter plot for Johansson, Krug and Mertins"""

    jo_pro_pro_scatter_plot = figure(title='NDUFS3 vs. NDUFS2', x_axis_label='protein z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")
    kr_pro_pro_scatter_plot = figure(title='NDUFS3 vs. NDUFS2', x_axis_label='protein z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")
    me_pro_pro_scatter_plot = figure(title='NDUFS3 vs. NDUFS2', x_axis_label='protein z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")

    return jo_pro_pro_scatter_plot, kr_pro_pro_scatter_plot, me_pro_pro_scatter_plot

def RNA_RNA_Scatter_Plot():
    """RNA-RNA for gene 1 and gene 2 scatter plot for Johansson, Krug and Mertins"""

    jo_rna_rna_scatter_plot = figure(title='NDUFS3 vs. NDUFS2', x_axis_label='mRNA z-score', y_axis_label='mRNA z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")
    kr_rna_rna_scatter_plot = figure(title='NDUFS3 vs. NDUFS2', x_axis_label='mRNA z-score', y_axis_label='mRNA z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")
    me_rna_rna_scatter_plot = figure(title='NDUFS3 vs. NDUFS2', x_axis_label='mRNA z-score', y_axis_label='mRNA z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")

    return jo_rna_rna_scatter_plot, kr_rna_rna_scatter_plot, me_rna_rna_scatter_plot

def Johansson_RNA_Pro_Scatter_Plot():
    """ mRNA-Protein scatter plot figures initiation"""

    jo_plot_mRNA_prot1 = figure(title='NDUFS2', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")
    jo_plot_mRNA_prot2 = figure(title='NDUFS3', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")
    jo_plot_mRNA_prot3 = figure(title='NDUFS7', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")
    jo_plot_mRNA_prot4 = figure(title='NDUFS8', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")

    return jo_plot_mRNA_prot1, jo_plot_mRNA_prot2, jo_plot_mRNA_prot3, jo_plot_mRNA_prot4


def Krug_RNA_Pro_Scatter_Plot():
    """ mRNA-Protein scatter plot figures initiation"""

    kr_plot_mRNA_prot1 = figure(title='NDUFS2', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")
    kr_plot_mRNA_prot2 = figure(title='NDUFS3', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")
    kr_plot_mRNA_prot3 = figure(title='NDUFS7', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")
    kr_plot_mRNA_prot4 = figure(title='NDUFS8', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")

    return kr_plot_mRNA_prot1, kr_plot_mRNA_prot2, kr_plot_mRNA_prot3, kr_plot_mRNA_prot4


def Mertins_RNA_Pro_Scatter_Plot():
    """ mRNA-Protein scatter plot figures initiation"""

    me_plot_mRNA_prot1 = figure(title='NDUFS2', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")

    me_plot_mRNA_prot2 = figure(title='NDUFS3', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")

    me_plot_mRNA_prot3 = figure(title='NDUFS7', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")

    me_plot_mRNA_prot4 = figure(title='NDUFS8', x_axis_label='mRNA z-score', y_axis_label='protein z-score',
                                width=SCATTER_PLOT_WIDTH, height=SCATTER_PLOT_HEIGHT, output_backend="webgl")

    return me_plot_mRNA_prot1, me_plot_mRNA_prot2, me_plot_mRNA_prot3, me_plot_mRNA_prot4