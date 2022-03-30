from bokeh.models import FactorRange
from bokeh.plotting import figure

from Plot_Line_Scatter_ColumnDataSource import Johansson_CDS, Krug_CDS, Mertins_CDS

# function call
[johansson_cds, johansson_subtype_tumor_tuple, jo_unique_gene_list] = Johansson_CDS()
[krug_cds, krug_subtype_tumor_tuple, kr_unique_gene_list] = Krug_CDS()
[mertins_cds, mertins_subtype_tumor_tuple, me_unique_gene_list] = Mertins_CDS()

# dictionary key name
subtype_tuple = 'subtype, tumor'
subtype = 'subtype'
protein_data = 'protein_data'
mRNA_data = 'mRNA_data'
gene = 'gene'

# figure size and other constants
LINE_PLOT_WIDTH = 1000
LINE_PLOT_HEIGHT = 200
GENE_NUMBER = ['Gene 1', 'Gene 2', 'Gene 3', 'Gene 4']
INITIAL_GENE = ['NDUFS2', 'NDUFS3', 'NDUFS7', 'blank 4']
GENE_COLORS = ['red', 'blue', 'green', 'orange']

def Johansson_Line_Plot():
    # Initializing the figure object
    jo_plot_p = figure(x_range=FactorRange(*johansson_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='protein z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)
    jo_plot_m = figure(x_range=FactorRange(*johansson_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='mRNA z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)

    for j in range(4):
        jo_plot_p.line(subtype_tuple, protein_data, source=johansson_cds[j], color=GENE_COLORS[j],
                       legend_label=GENE_NUMBER[j])
        jo_plot_p.circle(subtype_tuple, protein_data, source=johansson_cds[j], color=GENE_COLORS[j], size=4,
                         legend_label=GENE_NUMBER[j])
        jo_plot_m.line(subtype_tuple, mRNA_data, source=johansson_cds[j], color=GENE_COLORS[j],
                       legend_label=GENE_NUMBER[j])
        jo_plot_m.circle(subtype_tuple, mRNA_data, source=johansson_cds[j], color=GENE_COLORS[j], size=4,
                         legend_label=GENE_NUMBER[j])

    return jo_plot_p, jo_plot_m


def Krug_Line_Plot():
    kr_plot_p = figure(x_range=FactorRange(*krug_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='protein z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)
    kr_plot_m = figure(x_range=FactorRange(*krug_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='mRNA z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)

    for j in range(4):
        kr_plot_p.line(subtype_tuple, protein_data, source=krug_cds[j], color=GENE_COLORS[j],
                       legend_label=GENE_NUMBER[j])
        kr_plot_p.circle(subtype_tuple, protein_data, source=krug_cds[j], color=GENE_COLORS[j], size=4,
                         legend_label=GENE_NUMBER[j])
        kr_plot_m.line(subtype_tuple, mRNA_data, source=krug_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
        kr_plot_m.circle(subtype_tuple, mRNA_data, source=krug_cds[j], color=GENE_COLORS[j], size=4,
                         legend_label=GENE_NUMBER[j])

    return kr_plot_p, kr_plot_m


def Mertins_Line_Plot():
    me_plot_p = figure(x_range=FactorRange(*mertins_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='protein z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)
    me_plot_m = figure(x_range=FactorRange(*mertins_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='mRNA z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)

    for j in range(4):
        me_plot_p.line(subtype_tuple, protein_data, source=mertins_cds[j], color=GENE_COLORS[j],
                       legend_label=GENE_NUMBER[j])
        me_plot_p.circle(subtype_tuple, protein_data, source=mertins_cds[j], color=GENE_COLORS[j], size=4,
                         legend_label=GENE_NUMBER[j])
        me_plot_m.line(subtype_tuple, mRNA_data, source=mertins_cds[j], color=GENE_COLORS[j],
                       legend_label=GENE_NUMBER[j])
        me_plot_m.circle(subtype_tuple, mRNA_data, source=mertins_cds[j], color=GENE_COLORS[j], size=4,
                         legend_label=GENE_NUMBER[j])

    return me_plot_p, me_plot_m