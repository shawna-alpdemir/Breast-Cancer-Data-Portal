# function that style all plots
from bokeh.models import DataRange1d, ColorBar, LinearColorMapper, BasicTicker


def StylePlots(plot, PlotID):
    """ Aesthetic settings"""
    if PlotID == 'Correlation':
        plot.toolbar.autohide = True
        plot.xgrid.grid_line_color = None
        plot.ygrid.grid_line_color = None
        plot.x_range.range_padding = 0.05
        plot.xaxis.major_tick_line_color = None
        plot.xaxis.major_label_text_alpha = 0
        plot.yaxis.minor_tick_line_color = None
        plot.yaxis.ticker.min_interval = 2

        plot.legend.location = 'center'
        plot.legend.click_policy = 'hide'
        plot.add_layout(plot.legend[0], 'right')

    elif PlotID == 'mRNA-Prot':
        plot.title.align = 'center'
        plot.toolbar.autohide = True

        # plot.legend.label_text_font_size = '8pt'
        # plot.legend.location = "top_center"
        # plot.legend.orientation = "horizontal"
        # plot.legend.background_fill_alpha = 0
        # plot.legend.border_line_alpha = 0
        # plot.legend.spacing = 5
        # plot.legend.margin = -10
        # plot.legend.label_standoff = 3
        # plot.legend.glyph_width = 10
        # plot.add_layout(plot.legend[0], 'above')

        plot.xaxis.major_label_text_font_size = '8pt'
        plot.xaxis.minor_tick_line_color = None
        plot.xaxis.ticker.min_interval = 1
        plot.yaxis.minor_tick_line_color = None
        plot.yaxis.ticker.min_interval = 1
        plot.xgrid.grid_line_color = None
        plot.ygrid.grid_line_color = None


    elif PlotID == 'Subtype':
        plot.toolbar.autohide = True

        # plot.legend.label_text_font_size = '8pt'
        # plot.legend.location = "top_center"
        # plot.legend.orientation = "horizontal"
        # plot.legend.background_fill_alpha = 0
        # plot.legend.border_line_alpha = 0
        # plot.legend.spacing = 5
        # plot.legend.margin = -5
        # plot.legend.label_standoff = 3
        # plot.legend.glyph_width = 10
        # plot.legend.label_standoff = 0
        # plot.add_layout(plot.legend[0], 'above')

        plot.x_range.range_padding = 0.05
        plot.xaxis.major_tick_line_color = None
        plot.xaxis.major_label_text_alpha = 0
        #plot.xaxis.separator_line_color = None
        plot.yaxis.minor_tick_line_color = None
        plot.yaxis.ticker.min_interval = 1
        plot.xgrid.grid_line_color = None
        plot.ygrid.grid_line_color = None

    elif PlotID == 'Heatmap':
        plot.toolbar.autohide = True
        plot.grid.grid_line_color = None
        plot.axis.axis_line_color = None
        plot.axis.major_tick_line_color = None
        plot.xaxis.major_label_text_font_size = "8px"
        plot.yaxis.major_label_text_font_size = "10px"
        plot.axis.major_label_standoff = 0
        plot.xaxis.major_label_orientation = 3.14 / 3

        # set up color map, color using Sunset diverging colour scheme from https://personal.sron.nl/~pault/. this palette is related to RdYlBu scheme
        colors = ['#364B9A', '#4A7BB7', '#6EA6CD', '#98CAE1', '#C2E4EF', '#EAECCC', '#FEDA8B', '#FDB366', '#F67E4B',
                  '#DD3D2D', '#A50026']
        mapper = LinearColorMapper(palette=colors, low=-5, high=5)

        color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="10px",
                             ticker=BasicTicker(desired_num_ticks=len(colors)),
                             label_standoff=6, border_line_color=None)
        plot.add_layout(color_bar, 'right')

    return plot
