from bokeh.models import Whisker, FactorRange
from bokeh.plotting import figure
from bokeh.transform import factor_cmap
from functools import lru_cache # a function that helps in reducing the execution time of the function by using memoization technique

#testing
# import Import_data as im
# import Subtype_Average_DF as su
# import Subtype_Average_Plot_Source as ps
# import Style_Plot as style

# Global variables
SUBTYPE_PLOT_WIDTH = 500
SUBTYPE_PLOT_HEIGHT = 250

########################################################################################################################
################################################### Johansson ##########################################################

def jo_Subtype_Plot(jo_source_dict_subtype_protein, jo_source_dict_subtype_mrna, jo_subtype_colors, jo_subtypes):
    """ Function that carries bokeh plot configurations"""
    # Initialize the figure object
    jo_subtype_plot_protein = figure(x_range=FactorRange(*jo_source_dict_subtype_protein.data['x']),
                                     title='', x_axis_label='', y_axis_label='protein z-score', plot_height=SUBTYPE_PLOT_HEIGHT,
                                     plot_width=SUBTYPE_PLOT_WIDTH)

    jo_subtype_plot_mrna = figure(x_range=FactorRange(*jo_source_dict_subtype_mrna.data['x']),
                                  title='', x_axis_label='', y_axis_label='mRNA z-score', plot_height=SUBTYPE_PLOT_HEIGHT,
                                  plot_width=SUBTYPE_PLOT_WIDTH)

    # make protein plot
    jo_subtype_plot_protein.circle('x', 'y', source=jo_source_dict_subtype_protein,
                                   fill_color=factor_cmap('x', palette=jo_subtype_colors, factors=jo_subtypes, start=1,
                                                          end=2),
                                   line_color=factor_cmap('x', palette=jo_subtype_colors, factors=jo_subtypes, start=1,
                                                          end=2), legend_group='legend_group')
    jo_subtype_plot_protein.circle('x', 'upper', source=jo_source_dict_subtype_protein,
                                   size=0)  # invisible, here only to set initial visible range correctly
    jo_subtype_plot_protein.circle('x', 'lower', source=jo_source_dict_subtype_protein,
                                   size=0)  # invisible, here only to set initial visible range correctly

    # make mRNA plot
    jo_subtype_plot_mrna.circle('x', 'y', source=jo_source_dict_subtype_mrna,
                                fill_color=factor_cmap('x', palette=jo_subtype_colors, factors=jo_subtypes, start=1,
                                                       end=2),
                                line_color=factor_cmap('x', palette=jo_subtype_colors, factors=jo_subtypes, start=1,
                                                       end=2), legend_group='legend_group')
    jo_subtype_plot_mrna.circle('x', 'upper', source=jo_source_dict_subtype_mrna,
                                size=0)  # invisible, here only to set initial visible range correctly
    jo_subtype_plot_mrna.circle('x', 'lower', source=jo_source_dict_subtype_mrna,
                                size=0)  # invisible, here only to set initial visible range correctly

    # Add error bars
    # protein plot
    Whisker_protein = Whisker(source=jo_source_dict_subtype_protein, base="x", upper="upper", lower="lower",
                              line_color=factor_cmap('x', palette=jo_subtype_colors, factors=jo_subtypes, start=1,
                                                     end=2))
    Whisker_protein.upper_head.line_color = factor_cmap('x', palette=jo_subtype_colors, factors=jo_subtypes, start=1,
                                                        end=2)
    Whisker_protein.lower_head.line_color = factor_cmap('x', palette=jo_subtype_colors, factors=jo_subtypes, start=1,
                                                        end=2)
    jo_subtype_plot_protein.add_layout(Whisker_protein)

    # mRNA plot
    Whisker_mrna = Whisker(source=jo_source_dict_subtype_mrna, base="x", upper="upper", lower="lower",
                           line_color=factor_cmap('x', palette=jo_subtype_colors, factors=jo_subtypes, start=1, end=2))
    Whisker_mrna.upper_head.line_color = factor_cmap('x', palette=jo_subtype_colors, factors=jo_subtypes, start=1,
                                                     end=2)
    Whisker_mrna.lower_head.line_color = factor_cmap('x', palette=jo_subtype_colors, factors=jo_subtypes, start=1,
                                                     end=2)
    jo_subtype_plot_mrna.add_layout(Whisker_mrna)

    return jo_subtype_plot_protein, jo_subtype_plot_mrna


########################################################################################################################
##################################################### Krug #############################################################

def kr_Subtype_Plot(kr_source_dict_subtype_protein, kr_source_dict_subtype_mrna, kr_subtype_colors, kr_subtypes):
    """ Function that carries bokeh plot configurations"""
    # Initialize the figure object
    kr_subtype_plot_protein = figure(x_range=FactorRange(*kr_source_dict_subtype_protein.data['x']),
                                     title='', x_axis_label='', y_axis_label='protein z-score',
                                     plot_width=SUBTYPE_PLOT_WIDTH, plot_height=SUBTYPE_PLOT_HEIGHT)

    kr_subtype_plot_mrna = figure(x_range=FactorRange(*kr_source_dict_subtype_mrna.data['x']),
                                  title='', x_axis_label='', y_axis_label='mRNA z-score', plot_width=SUBTYPE_PLOT_WIDTH,
                                  plot_height=SUBTYPE_PLOT_HEIGHT)

    # make protein plots
    kr_subtype_plot_protein.circle('x', 'y', source=kr_source_dict_subtype_protein,
                                   fill_color=factor_cmap('x', palette=kr_subtype_colors, factors=kr_subtypes, start=1,
                                                          end=2),
                                   line_color=factor_cmap('x', palette=kr_subtype_colors, factors=kr_subtypes, start=1,
                                                          end=2), legend_group='legend_group')
    kr_subtype_plot_protein.circle('x', 'upper', source=kr_source_dict_subtype_protein,
                                   size=0)  # invisible, here only to set initial visible range correctly
    kr_subtype_plot_protein.circle('x', 'lower', source=kr_source_dict_subtype_protein,
                                   size=0)  # invisible, here only to set initial visible range correctly

    # make mRNA plots
    kr_subtype_plot_mrna.circle('x', 'y', source=kr_source_dict_subtype_mrna,
                                fill_color=factor_cmap('x', palette=kr_subtype_colors, factors=kr_subtypes, start=1,
                                                       end=2),
                                line_color=factor_cmap('x', palette=kr_subtype_colors, factors=kr_subtypes, start=1,
                                                       end=2), legend_group='legend_group')
    kr_subtype_plot_mrna.circle('x', 'upper', source=kr_source_dict_subtype_mrna,
                                size=0)  # invisible, here only to set initial visible range correctly
    kr_subtype_plot_mrna.circle('x', 'lower', source=kr_source_dict_subtype_mrna,
                                size=0)  # invisible, here only to set initial visible range correctly

    # Add error bars
    # protein plot
    Whisker_protein = Whisker(source=kr_source_dict_subtype_protein, base="x", upper="upper", lower="lower",
                              line_color=factor_cmap('x', palette=kr_subtype_colors, factors=kr_subtypes, start=1,
                                                     end=2))
    Whisker_protein.upper_head.line_color = factor_cmap('x', palette=kr_subtype_colors, factors=kr_subtypes, start=1,
                                                        end=2)
    Whisker_protein.lower_head.line_color = factor_cmap('x', palette=kr_subtype_colors, factors=kr_subtypes, start=1,
                                                        end=2)
    kr_subtype_plot_protein.add_layout(Whisker_protein)

    # mRNA plot
    Whisker_mrna = Whisker(source=kr_source_dict_subtype_mrna, base="x", upper="upper", lower="lower",
                           line_color=factor_cmap('x', palette=kr_subtype_colors, factors=kr_subtypes, start=1, end=2))
    Whisker_mrna.upper_head.line_color = factor_cmap('x', palette=kr_subtype_colors, factors=kr_subtypes, start=1,
                                                     end=2)
    Whisker_mrna.lower_head.line_color = factor_cmap('x', palette=kr_subtype_colors, factors=kr_subtypes, start=1,
                                                     end=2)
    kr_subtype_plot_mrna.add_layout(Whisker_mrna)

    return kr_subtype_plot_protein, kr_subtype_plot_mrna


########################################################################################################################
##################################################### Mertins ##########################################################

def me_Subtype_Plot(me_source_dict_subtype_protein, me_source_dict_subtype_mrna, me_subtype_colors, me_subtypes):
    """ Function that carries bokeh plot configurations"""
    # Initialize the figure object
    me_subtype_plot_protein = figure(x_range=FactorRange(*me_source_dict_subtype_protein.data['x']),
                                     title='', x_axis_label='', y_axis_label='protein z-score',
                                     plot_width=SUBTYPE_PLOT_WIDTH, plot_height=SUBTYPE_PLOT_HEIGHT)

    me_subtype_plot_mrna = figure(x_range=FactorRange(*me_source_dict_subtype_mrna.data['x']),
                                  title='', x_axis_label='', y_axis_label='mRNA z-score', plot_width=SUBTYPE_PLOT_WIDTH,
                                  plot_height=SUBTYPE_PLOT_HEIGHT)

    # make protein plots
    me_subtype_plot_protein.circle('x', 'y', source=me_source_dict_subtype_protein,
                                   fill_color=factor_cmap('x', palette=me_subtype_colors, factors=me_subtypes, start=1,
                                                          end=2),
                                   line_color=factor_cmap('x', palette=me_subtype_colors, factors=me_subtypes, start=1,
                                                          end=2), legend_group='legend_group')
    me_subtype_plot_protein.circle('x', 'upper', source=me_source_dict_subtype_protein,
                                   size=0)  # invisible, here only to set initial visible range correctly
    me_subtype_plot_protein.circle('x', 'lower', source=me_source_dict_subtype_protein,
                                   size=0)  # invisible, here only to set initial visible range correctly

    # make mRNA plots
    me_subtype_plot_mrna.circle('x', 'y', source=me_source_dict_subtype_mrna,
                                fill_color=factor_cmap('x', palette=me_subtype_colors, factors=me_subtypes, start=1,
                                                       end=2),
                                line_color=factor_cmap('x', palette=me_subtype_colors, factors=me_subtypes, start=1,
                                                       end=2), legend_group='legend_group')
    me_subtype_plot_mrna.circle('x', 'upper', source=me_source_dict_subtype_mrna,
                                size=0)  # invisible, here only to set initial visible range correctly
    me_subtype_plot_mrna.circle('x', 'lower', source=me_source_dict_subtype_mrna,
                                size=0)  # invisible, here only to set initial visible range correctly

    # Add error bars
    # protein plot
    Whisker_protein = Whisker(source=me_source_dict_subtype_protein, base="x", upper="upper", lower="lower",
                              line_color=factor_cmap('x', palette=me_subtype_colors, factors=me_subtypes, start=1,
                                                     end=2))
    Whisker_protein.upper_head.line_color = factor_cmap('x', palette=me_subtype_colors, factors=me_subtypes, start=1,
                                                        end=2)
    Whisker_protein.lower_head.line_color = factor_cmap('x', palette=me_subtype_colors, factors=me_subtypes, start=1,
                                                        end=2)
    me_subtype_plot_protein.add_layout(Whisker_protein)

    # mRNA plot
    Whisker_mrna = Whisker(source=me_source_dict_subtype_mrna, base="x", upper="upper", lower="lower",
                           line_color=factor_cmap('x', palette=me_subtype_colors, factors=me_subtypes, start=1, end=2))
    Whisker_mrna.upper_head.line_color = factor_cmap('x', palette=me_subtype_colors, factors=me_subtypes, start=1,
                                                     end=2)
    Whisker_mrna.lower_head.line_color = factor_cmap('x', palette=me_subtype_colors, factors=me_subtypes, start=1,
                                                     end=2)
    me_subtype_plot_mrna.add_layout(Whisker_mrna)

    return me_subtype_plot_protein, me_subtype_plot_mrna

# testing
# # return jo_df_protein, jo_df_mrna, jo_genes, jo_subtypes, jo_subtype_colors
# a = im.jo_ImportData()
#
# # return jo_protein_mean_by_subtype, jo_mrna_mean_by_subtype, jo_protein_sem_by_subtype, jo_mrna_sem_by_subtype
# b = su.jo_Subtype_Avg_SEM_DFs(a[0],a[1],a[3])
#
# # return jo_source_dict_subtype_protein, jo_source_dict_subtype_mrna, jo_subtypes
# c = ps.jo_Subtype_Avg_Plot_Source(b[0],b[1],b[2],b[3],a[3], ['ESR1','PGR','ERBB2','MKI67'])
#
# # return jo_subtype_plot_protein, jo_subtype_plot_mrna
# d = jo_Subtype_Plot(c[0], c[1], a[4], a[3])
#
# # show the test plot
# show(style.StylePlots(d[1], 'Subtype'))
