# function that style all plots
from bokeh.models import DataRange1d


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

    return plot
