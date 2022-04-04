from bokeh.models import Legend, DataRange1d


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
        plot.legend.location = 'center'
        plot.add_layout(plot.legend[0], 'right')
        plot.legend.click_policy = 'hide'

    elif PlotID == 'mRNA-Prot':
        plot.xaxis.major_label_text_font_size = '8pt'
        plot.title.align = 'center'
        plot.toolbar.autohide = True
        plot.add_layout(plot.legend[0], 'right')
        plot.xaxis.minor_tick_line_color = None
        plot.yaxis.minor_tick_line_color = None
        plot.xgrid.grid_line_color = None
        plot.ygrid.grid_line_color = None

    elif PlotID == 'Subtype':
        plot.legend.label_standoff = 0
        plot.toolbar.autohide = True
        plot.add_layout(plot.legend[0], 'right')
        plot.x_range.range_padding = 0.05
        plot.xaxis.major_tick_line_color = None
        plot.xaxis.major_label_text_alpha = 0
        plot.yaxis.minor_tick_line_color = None
        #plot.xaxis.separator_line_color = None
        plot.xgrid.grid_line_color = None
        plot.ygrid.grid_line_color = None

    return plot
