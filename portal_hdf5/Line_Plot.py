# function that initiates the line plot figure

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
LINE_PLOT_HEIGHT = 125

def Johansson_Line_Plot():
    """Initializing the line plot"""
    jo_plot_p = figure(x_range=FactorRange(*johansson_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='protein z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)
    jo_plot_m = figure(x_range=FactorRange(*johansson_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='mRNA z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)

    return jo_plot_p, jo_plot_m


def Krug_Line_Plot():
    """Initializing the line plot"""
    kr_plot_p = figure(x_range=FactorRange(*krug_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='protein z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)
    kr_plot_m = figure(x_range=FactorRange(*krug_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='mRNA z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)

    return kr_plot_p, kr_plot_m


def Mertins_Line_Plot():
    """Initializing the line plot"""
    me_plot_p = figure(x_range=FactorRange(*mertins_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='protein z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)
    me_plot_m = figure(x_range=FactorRange(*mertins_subtype_tumor_tuple), title='', x_axis_label='',
                       y_axis_label='mRNA z-score', plot_width=LINE_PLOT_WIDTH, plot_height=LINE_PLOT_HEIGHT)

    return me_plot_p, me_plot_m