# function that returns download button for gene entry, download button for correlation table is in main.py

# import
from os.path import join, dirname
from bokeh.models import CustomJS, Button
from Plot_Line_Scatter_ColumnDataSource import Johansson_CDS, Krug_CDS, Mertins_CDS

# function calls
[johansson_cds, johansson_subtype_tumor_tuple, jo_unique_gene_list] = Johansson_CDS()
[krug_cds, krug_subtype_tumor_tuple, kr_unique_gene_list] = Krug_CDS()
[mertins_cds, mertins_subtype_tumor_tuple, me_unique_gene_list] = Mertins_CDS()

def Download_Buttons():
    button1 = Button(label="Download Gene 1 Data", button_type="success", width = 150)
    button1.js_on_event("button_click", CustomJS(args=dict(source=johansson_cds[0]), code=open(join(dirname(__file__),
                                                                                                    "Download_Javascript/jo_download.js")).read()))
    button1.js_on_event("button_click", CustomJS(args=dict(source=krug_cds[0]), code=open(join(dirname(__file__),
                                                                                               "Download_Javascript/kr_download.js")).read()))
    button1.js_on_event("button_click", CustomJS(args=dict(source=mertins_cds[0]), code=open(join(dirname(__file__),
                                                                                                  "Download_Javascript/me_download.js")).read()))
    ####
    button2 = Button(label="Download Gene 2 Data", button_type="success", width = 150)
    button2.js_on_event("button_click", CustomJS(args=dict(source=johansson_cds[1]), code=open(join(dirname(__file__),
                                                                                                    "Download_Javascript/jo_download.js")).read()))
    button2.js_on_event("button_click", CustomJS(args=dict(source=krug_cds[1]), code=open(join(dirname(__file__),
                                                                                               "Download_Javascript/kr_download.js")).read()))
    button2.js_on_event("button_click", CustomJS(args=dict(source=mertins_cds[1]), code=open(join(dirname(__file__),
                                                                                                  "Download_Javascript/me_download.js")).read()))
    ####
    button3 = Button(label="Download Gene 3 Data", button_type="success", width = 150)
    button3.js_on_event("button_click", CustomJS(args=dict(source=johansson_cds[2]), code=open(join(dirname(__file__),
                                                                                                    "Download_Javascript/jo_download.js")).read()))
    button3.js_on_event("button_click", CustomJS(args=dict(source=krug_cds[2]), code=open(join(dirname(__file__),
                                                                                               "Download_Javascript/kr_download.js")).read()))
    button3.js_on_event("button_click", CustomJS(args=dict(source=mertins_cds[2]), code=open(join(dirname(__file__),
                                                                                                  "Download_Javascript/me_download.js")).read()))
    ####
    button4 = Button(label="Download Gene 4 Data", button_type="success", width = 150)
    button4.js_on_event("button_click", CustomJS(args=dict(source=johansson_cds[3]), code=open(join(dirname(__file__),
                                                                                                    "Download_Javascript/jo_download.js")).read()))
    button4.js_on_event("button_click", CustomJS(args=dict(source=krug_cds[3]), code=open(join(dirname(__file__),
                                                                                               "Download_Javascript/kr_download.js")).read()))
    button4.js_on_event("button_click", CustomJS(args=dict(source=mertins_cds[3]), code=open(join(dirname(__file__),
                                                                                                  "Download_Javascript/me_download.js")).read()))

    return button1, button2, button3, button4