# function that take subtype plot CDS to make subtype plots

#from pdb import set_trace
from bokeh.models import Whisker, FactorRange
from bokeh.plotting import figure
from bokeh.transform import factor_cmap
#from Subtype_Average_Plot_CDS import Johansson_Subtype_Avg_Plot_CDS, Krug_Subtype_Avg_Plot_CDS, Mertins_Subtype_Avg_Plot_CDS
from Subtype_Plot_CDS_new import Johansson_Subtype_Avg_Plot_CDS, Krug_Subtype_Avg_Plot_CDS, Mertins_Subtype_Avg_Plot_CDS

# constants
SUBTYPE_PLOT_WIDTH = 410
SUBTYPE_PLOT_HEIGHT = 250
five_subtype_colors = ['#E31A1C', '#FB9A99', '#1F78B4', '#A6CEE3','#33A02C']
five_subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']

# function call
[jo_subtype_protein_CDS, jo_subtype_mRNA_CDS] = Johansson_Subtype_Avg_Plot_CDS()
[kr_subtype_protein_CDS, kr_subtype_mRNA_CDS] = Krug_Subtype_Avg_Plot_CDS()
[me_subtype_protein_CDS, me_subtype_mRNA_CDS] = Mertins_Subtype_Avg_Plot_CDS()

########################################################################################################################
################################################### Johansson ##########################################################

def Johansson_Subtype_Plot(jo_subtype_protein_CDS, jo_subtype_mRNA_CDS):
    """ Function that carries bokeh plot configurations"""
    # Initialize the figure object
    jo_protein_subtype_plot = figure(x_range=FactorRange(*jo_subtype_protein_CDS.data["x"]),
                                     title='', x_axis_label='', y_axis_label='protein z-score',
                                     height=SUBTYPE_PLOT_HEIGHT, width=SUBTYPE_PLOT_WIDTH, output_backend="webgl")

    jo_mRNA_subtype_plot = figure(x_range=FactorRange(*jo_subtype_mRNA_CDS.data['x']),
                                  title='', x_axis_label='', y_axis_label='mRNA z-score',
                                  height=SUBTYPE_PLOT_HEIGHT, width=SUBTYPE_PLOT_WIDTH, output_backend="webgl")

    # make protein plot
    jo_protein_subtype_plot.scatter('x', 'y', source=jo_subtype_protein_CDS,
                                   fill_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2),
                                   line_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2),
                                   )
    jo_protein_subtype_plot.scatter('x', 'upper', source=jo_subtype_protein_CDS, size=0)  # invisible, here only to set initial visible range correctly
    jo_protein_subtype_plot.scatter('x', 'lower', source=jo_subtype_protein_CDS, size=0)  # invisible, here only to set initial visible range correctly

    # make mRNA plot
    jo_mRNA_subtype_plot.scatter('x', 'y', source=jo_subtype_mRNA_CDS,
                                fill_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2),
                                line_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2),
                                )
    jo_mRNA_subtype_plot.scatter('x', 'upper', source=jo_subtype_mRNA_CDS, size=0)  # invisible, here only to set initial visible range correctly
    jo_mRNA_subtype_plot.scatter('x', 'lower', source=jo_subtype_mRNA_CDS, size=0)  # invisible, here only to set initial visible range correctly

    # Add error bars
    # protein plot
    Whisker_protein = Whisker(source=jo_subtype_protein_CDS, base="x", upper="upper", lower="lower",
                              line_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2))
    Whisker_protein.upper_head.line_color = factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2)
    Whisker_protein.lower_head.line_color = factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2)
    jo_protein_subtype_plot.add_layout(Whisker_protein)

    # mRNA plot
    Whisker_mrna = Whisker(source=jo_subtype_mRNA_CDS, base="x", upper="upper", lower="lower",
                           line_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2))
    Whisker_mrna.upper_head.line_color = factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2)
    Whisker_mrna.lower_head.line_color = factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2)
    jo_mRNA_subtype_plot.add_layout(Whisker_mrna)

    return jo_protein_subtype_plot, jo_mRNA_subtype_plot


########################################################################################################################
##################################################### Krug #############################################################

def Krug_Subtype_Plot(kr_subtype_protein_CDS, kr_subtype_mRNA_CDS):
    """ Function that carries bokeh plot configurations"""
    # Initialize the figure object
    kr_protein_subtype_plot = figure(x_range=FactorRange(*kr_subtype_protein_CDS.data["x"]),
                                     title='', x_axis_label='', y_axis_label='protein z-score',
                                     height=SUBTYPE_PLOT_HEIGHT, width=SUBTYPE_PLOT_WIDTH, output_backend="webgl")

    kr_mRNA_subtype_plot = figure(x_range=FactorRange(*kr_subtype_mRNA_CDS.data['x']),
                                  title='', x_axis_label='', y_axis_label='mRNA z-score',
                                  height=SUBTYPE_PLOT_HEIGHT, width=SUBTYPE_PLOT_WIDTH, output_backend="webgl")

    # make protein plot
    kr_protein_subtype_plot.scatter('x', 'y', source=kr_subtype_protein_CDS,
                                   fill_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2),
                                   line_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2),
                                   )
    kr_protein_subtype_plot.scatter('x', 'upper', source=kr_subtype_protein_CDS, size=0)  # invisible, here only to set initial visible range correctly
    kr_protein_subtype_plot.scatter('x', 'lower', source=kr_subtype_protein_CDS, size=0)  # invisible, here only to set initial visible range correctly

    # make mRNA plot
    kr_mRNA_subtype_plot.scatter('x', 'y', source=kr_subtype_mRNA_CDS,
                                fill_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2),
                                line_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2),
                                )
    kr_mRNA_subtype_plot.scatter('x', 'upper', source=kr_subtype_mRNA_CDS, size=0)  # invisible, here only to set initial visible range correctly
    kr_mRNA_subtype_plot.scatter('x', 'lower', source=kr_subtype_mRNA_CDS, size=0)  # invisible, here only to set initial visible range correctly

    # Add error bars
    # protein plot
    Whisker_protein = Whisker(source=kr_subtype_protein_CDS, base="x", upper="upper", lower="lower",
                              line_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2))
    Whisker_protein.upper_head.line_color = factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2)
    Whisker_protein.lower_head.line_color = factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2)
    kr_protein_subtype_plot.add_layout(Whisker_protein)

    # mRNA plot
    Whisker_mrna = Whisker(source=kr_subtype_mRNA_CDS, base="x", upper="upper", lower="lower",
                           line_color=factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2))
    Whisker_mrna.upper_head.line_color = factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2)
    Whisker_mrna.lower_head.line_color = factor_cmap('x', palette=five_subtype_colors, factors=five_subtypes, start=1, end=2)
    kr_mRNA_subtype_plot.add_layout(Whisker_mrna)

    return kr_protein_subtype_plot, kr_mRNA_subtype_plot


########################################################################################################################
##################################################### Mertins ##########################################################

def Mertins_Subtype_Plot(me_subtype_protein_CDS, me_subtype_mRNA_CDS):
    """ Function that carries bokeh plot configurations"""
    # Initialize the figure object
    me_protein_subtype_plot = figure(x_range=FactorRange(*me_subtype_protein_CDS.data["x"]),
                                     title='', x_axis_label='', y_axis_label='protein z-score',
                                     height=SUBTYPE_PLOT_HEIGHT, width=SUBTYPE_PLOT_WIDTH, output_backend="webgl")

    me_mRNA_subtype_plot = figure(x_range=FactorRange(*me_subtype_mRNA_CDS.data['x']),
                                  title='', x_axis_label='', y_axis_label='mRNA z-score',
                                  height=SUBTYPE_PLOT_HEIGHT, width=SUBTYPE_PLOT_WIDTH, output_backend="webgl")

    # make protein plot
    me_protein_subtype_plot.scatter('x', 'y', source=me_subtype_protein_CDS,
                                   fill_color=factor_cmap('x', palette=five_subtype_colors[0:4], factors=five_subtypes[0:4], start=1, end=2),
                                   line_color=factor_cmap('x', palette=five_subtype_colors[0:4], factors=five_subtypes[0:4], start=1, end=2),
                                   )
    me_protein_subtype_plot.scatter('x', 'upper', source=me_subtype_protein_CDS, size=0)  # invisible, here only to set initial visible range correctly
    me_protein_subtype_plot.scatter('x', 'lower', source=me_subtype_protein_CDS, size=0)  # invisible, here only to set initial visible range correctly

    # make mRNA plot
    me_mRNA_subtype_plot.scatter('x', 'y', source=me_subtype_mRNA_CDS,
                                fill_color=factor_cmap('x', palette=five_subtype_colors[0:4], factors=five_subtypes[0:4], start=1, end=2),
                                line_color=factor_cmap('x', palette=five_subtype_colors[0:4], factors=five_subtypes[0:4], start=1, end=2),
                                )
    me_mRNA_subtype_plot.scatter('x', 'upper', source=me_subtype_mRNA_CDS, size=0)  # invisible, here only to set initial visible range correctly
    me_mRNA_subtype_plot.scatter('x', 'lower', source=me_subtype_mRNA_CDS, size=0)  # invisible, here only to set initial visible range correctly

    # Add error bars
    # protein plot
    Whisker_protein = Whisker(source=me_subtype_protein_CDS, base="x", upper="upper", lower="lower",
                              line_color=factor_cmap('x', palette=five_subtype_colors[0:4], factors=five_subtypes[0:4], start=1, end=2))
    Whisker_protein.upper_head.line_color = factor_cmap('x', palette=five_subtype_colors[0:4], factors=five_subtypes[0:4], start=1, end=2)
    Whisker_protein.lower_head.line_color = factor_cmap('x', palette=five_subtype_colors[0:4], factors=five_subtypes[0:4], start=1, end=2)
    me_protein_subtype_plot.add_layout(Whisker_protein)

    # mRNA plot
    Whisker_mrna = Whisker(source=me_subtype_mRNA_CDS, base="x", upper="upper", lower="lower",
                           line_color=factor_cmap('x', palette=five_subtype_colors[0:4], factors=five_subtypes[0:4], start=1, end=2))
    Whisker_mrna.upper_head.line_color = factor_cmap('x', palette=five_subtype_colors[0:4], factors=five_subtypes[0:4], start=1, end=2)
    Whisker_mrna.lower_head.line_color = factor_cmap('x', palette=five_subtype_colors[0:4], factors=five_subtypes[0:4], start=1, end=2)
    me_mRNA_subtype_plot.add_layout(Whisker_mrna)

    return me_protein_subtype_plot, me_mRNA_subtype_plot