# --------------------------------------------------------------------
# Import packages and scripts
# --------------------------------------------------------------------
import h5py
import numpy as np
import time
from os.path import join, dirname
from bokeh.io import curdoc
from bokeh.layouts import layout, column, row
from bokeh.models import Spacer, Tabs, TabPanel, AutocompleteInput, ColumnDataSource, TableColumn, DataTable, Button, \
    Select, Div, CustomJS, ScientificFormatter
from bokeh.transform import factor_cmap
from bokeh.embed import server_document
from Correlation_Table import Get_Protein_Correlation_Table, Get_mRNA_Correlation_Table
from Gene_List import Gene_List
from Import_Files import Import_HDF5, Import_Static_Correlation_Table, Import_Subtype_DF
from Plot_All_Styling import StylePlots
from Plot_Line_Scatter_ColumnDataSource import Johansson_CDS, Krug_CDS, Mertins_CDS, Johansson_Scatter_CDS, Krug_Scatter_CDS, Mertins_Scatter_CDS
from Line_Plot import Johansson_Line_Plot, Krug_Line_Plot, Mertins_Line_Plot
from Scatter_Plot import Pro_Pro_Scatter_Plot, RNA_RNA_Scatter_Plot, Johansson_RNA_Pro_Scatter_Plot, Krug_RNA_Pro_Scatter_Plot, Mertins_RNA_Pro_Scatter_Plot
from Subtype_Plot_CDS_new import Johansson_Subtype_Avg_Plot_CDS, Krug_Subtype_Avg_Plot_CDS, Mertins_Subtype_Avg_Plot_CDS
from Subtype_Plot import Johansson_Subtype_Plot, Krug_Subtype_Plot, Mertins_Subtype_Plot

# --------------------------------------------------------------------
#Row 0: Functions Def
# --------------------------------------------------------------------

# constants
TICKER = {}
TICKER_INDEX = None
SUBTYPE_PLOT_GENE_LIST = ['ESR1', 'PGR', 'ERBB2', 'MKI67']
TEXTBOX_GENE_LIST = ['NDUFS2', 'NDUFS3', 'NDUFS7', 'NDUFS8']

# dictionary key name
subtype_tuple = 'subtype, tumor'
subtype = 'subtype'
protein_data = 'protein_data'
mRNA_data = 'mRNA_data'
gene = 'gene'

# function call
[JohanssonProteome, JohanssonTranscriptome, KrugProteome, KrugTranscriptome, MertinsProteome, MertinsTranscriptome] = Import_HDF5()
[jo_pro, jo_rna, kr_rna, kr_pro, me_rna, me_pro] = Import_Subtype_DF()

# --------------------------------------------------------------------
# Textbox value tracing
# --------------------------------------------------------------------
def Nix(val, lst):
    """function that prevents user input the same gene in other tickers"""
    return [x for x in lst if x!= val]

def Ticker1_Change(attrname, old, new):
    """ change function that triggers the update"""
    global TICKER_INDEX
    TICKER_INDEX = 0
    TICKER[1].completions = Nix(new, all_unique_genes)
    TICKER[2].completions = Nix(new, all_unique_genes)
    TICKER[3].completions = Nix(new, all_unique_genes)

def Status1_Change(attrname, old, new):
    gene_hide_text_div.text = f"Gene 1 has been changed to {new}"

def Ticker2_Change(attrname, old, new):
    global TICKER_INDEX
    TICKER_INDEX = 1
    TICKER[0].completions = Nix(new, all_unique_genes)
    TICKER[2].completions = Nix(new, all_unique_genes)
    TICKER[3].completions = Nix(new, all_unique_genes)

def Status2_Change(attrname, old, new):
    gene_hide_text_div.text = f"Gene 2 has been changed to {new}"

def Ticker3_Change(attrname, old, new):
    global TICKER_INDEX
    TICKER_INDEX = 2
    TICKER[0].completions = Nix(new, all_unique_genes)
    TICKER[1].completions = Nix(new, all_unique_genes)
    TICKER[3].completions = Nix(new, all_unique_genes)

def Status3_Change(attrname, old, new):
    gene_hide_text_div.text = f"Gene 3 has been changed to {new}"

def Ticker4_Change(attrname, old, new):
    global TICKER_INDEX
    TICKER_INDEX = 3
    TICKER[0].completions = Nix(new, all_unique_genes)
    TICKER[1].completions = Nix(new, all_unique_genes)
    TICKER[2].completions = Nix(new, all_unique_genes)

def Status4_Change(attrname, old, new):
    gene_hide_text_div.text = f"Gene 4 has been changed to {new}"

TICKER_FUNCTION_LIST = [Ticker1_Change, Ticker2_Change, Ticker3_Change, Ticker4_Change]
STATUS_FUNCTION_LIST = [Status1_Change, Status2_Change, Status3_Change, Status4_Change]

# --------------------------------------------------------------------
# Line and scatter plot (same CDS)
# --------------------------------------------------------------------
def Line_Plot_Update(attrname, old, new):
    """Update function for text entries call back targeting line plot"""

    with h5py.File(JohanssonProteome, "r") as a, h5py.File(JohanssonTranscriptome, "r") as b:
        protein_data_new_gene = np.array(a.get(new))
        mRNA_data_new_gene = np.array(b.get(new))

        if not np.all(protein_data_new_gene):
            protein_data_new_gene = np.repeat(np.nan, 45) # line plot takes NaN for None data

        if not np.all(mRNA_data_new_gene):
            mRNA_data_new_gene = np.repeat(np.nan, 45)

        johansson_cds[TICKER_INDEX].data = {subtype_tuple: johansson_subtype_tumor_tuple,
                 subtype: list(zip(*johansson_subtype_tumor_tuple))[0], # extract first element of the list of subtype_tumor_tuple, which is the subtype
                 protein_data: protein_data_new_gene,
                 mRNA_data: mRNA_data_new_gene,
                 gene: np.repeat(new, 45)}

    with h5py.File(KrugProteome, "r") as a, h5py.File(KrugTranscriptome, "r") as b:
        protein_data_new_gene = np.array(a.get(new))
        mRNA_data_new_gene = np.array(b.get(new))

        if not np.all(protein_data_new_gene):
            protein_data_new_gene = np.repeat(np.nan, 122)

        if not np.all(mRNA_data_new_gene):
            mRNA_data_new_gene = np.repeat(np.nan, 122)

        krug_cds[TICKER_INDEX].data = {subtype_tuple: krug_subtype_tumor_tuple,
                 subtype: list(zip(*krug_subtype_tumor_tuple))[0], # extract first element of the list of subtype_tumor_tuple
                 protein_data: protein_data_new_gene,
                 mRNA_data: mRNA_data_new_gene,
                 gene: np.repeat(new, 122)}

    with h5py.File(MertinsProteome, "r") as a, h5py.File(MertinsTranscriptome, "r") as b:
        protein_data_new_gene = np.array(a.get(new))
        mRNA_data_new_gene = np.array(b.get(new))

        if not np.all(protein_data_new_gene):
            protein_data_new_gene = np.repeat(np.nan, 77)

        if not np.all(mRNA_data_new_gene):
            mRNA_data_new_gene = np.repeat(np.nan, 77)

        mertins_cds[TICKER_INDEX].data = {subtype_tuple: mertins_subtype_tumor_tuple,
                 subtype: list(zip(*mertins_subtype_tumor_tuple))[0], # extract first element of the list of subtype_tumor_tuple
                 protein_data: protein_data_new_gene,
                 mRNA_data: mRNA_data_new_gene,
                 gene: np.repeat(new, 77)}

    print(f"Line plots updated...")

# --------------------------------------------------------------------
# Scatter plot update button
# --------------------------------------------------------------------
def Textbox_GeneList_Update(attrname, old, new):
    global TEXTBOX_GENE_LIST
    TEXTBOX_GENE_LIST[TICKER_INDEX] = new
    print(f'Textbox gene list updated: {TEXTBOX_GENE_LIST}...')

def Scatter_Plot_Update(event):
    """Update function for text entries call back targeting line plot"""

    for j, k in zip(TEXTBOX_GENE_LIST, range(4)):
        with h5py.File(JohanssonProteome, "r") as a, h5py.File(JohanssonTranscriptome, "r") as b:
            protein_data_new_gene = np.array(a.get(j))
            mRNA_data_new_gene = np.array(b.get(j))

            if not np.all(protein_data_new_gene):
                protein_data_new_gene = np.zeros(45)

            if not np.all(mRNA_data_new_gene):
                mRNA_data_new_gene = np.zeros(45)

            johansson_scatter_cds[k].data = {subtype: list(zip(*johansson_subtype_tumor_tuple))[0],
                                             # extract first element of the list of subtype_tumor_tuple, which is the subtype
                                             protein_data: protein_data_new_gene,
                                             mRNA_data: mRNA_data_new_gene}
            jo_plot_mRNA_prot[k].title.text = j

        with h5py.File(KrugProteome, "r") as a, h5py.File(KrugTranscriptome, "r") as b:
            protein_data_new_gene = np.array(a.get(j))
            mRNA_data_new_gene = np.array(b.get(j))

            if not np.all(protein_data_new_gene):
                protein_data_new_gene = np.zeros(122)

            if not np.all(mRNA_data_new_gene):
                mRNA_data_new_gene = np.zeros(122)

            krug_scatter_cds[k].data = {subtype: list(zip(*krug_subtype_tumor_tuple))[0],
                                        # extract first element of the list of subtype_tumor_tuple
                                        protein_data: protein_data_new_gene,
                                        mRNA_data: mRNA_data_new_gene}
            kr_plot_mRNA_prot[k].title.text = j

        with h5py.File(MertinsProteome, "r") as a, h5py.File(MertinsTranscriptome, "r") as b:
            protein_data_new_gene = np.array(a.get(j))
            mRNA_data_new_gene = np.array(b.get(j))

            if not np.all(protein_data_new_gene):
                protein_data_new_gene = np.zeros(77)

            if not np.all(mRNA_data_new_gene):
                mRNA_data_new_gene = np.zeros(77)

            mertins_scatter_cds[k].data = {subtype: list(zip(*mertins_subtype_tumor_tuple))[0],
                                                      # extract first element of the list of subtype_tumor_tuple
                                                      protein_data: protein_data_new_gene,
                                                      mRNA_data: mRNA_data_new_gene}
            me_plot_mRNA_prot[k].title.text = j

    print(f"Scatter plot updated...")

def Scatter_Select_Button_Update(event):
    """ Update function for select menu; it locates the value in textbox and return """
    option_list = ['Gene 1', 'Gene 2', 'Gene 3', 'Gene 4']
    x_index = option_list.index(x_axis_select.value) # get the index of which gene users selected
    y_index = option_list.index(y_axis_select.value)

    x_gene_name = TICKER[x_index].value # use the index to locate what gene is in the corresponding textbox
    y_gene_name = TICKER[y_index].value

    with h5py.File(JohanssonProteome, "r") as a, h5py.File(JohanssonTranscriptome, "r") as b:
        x_protein_data = np.array(a.get(x_gene_name))
        y_protein_data = np.array(a.get(y_gene_name))

        x_mrna_data = np.array(b.get(x_gene_name))
        y_mrna_data = np.array(b.get(y_gene_name))

        old_data_list = [x_protein_data,y_protein_data,x_mrna_data,y_mrna_data]
        new_data_list = []
        for data in old_data_list:
            if not np.all(data):
                data = np.zeros(45)
            new_data_list.append(data)

        johansson_cds[4].data = {subtype_tuple: johansson_subtype_tumor_tuple,
                                            subtype: list(zip(*johansson_subtype_tumor_tuple))[0],
                                            'x_protein_data': new_data_list[0], # index 0
                                            'y_protein_data': new_data_list[1]} # index 1
        johansson_cds[5].data = {subtype_tuple: johansson_subtype_tumor_tuple,
                                            subtype: list(zip(*johansson_subtype_tumor_tuple))[0],
                                            'x_mRNA_data': new_data_list[2], # index 2
                                            'y_mRNA_data': new_data_list[3]} # index 3

    with h5py.File(KrugProteome, "r") as a, h5py.File(KrugTranscriptome, "r") as b:
        x_protein_data = np.array(a.get(x_gene_name))
        y_protein_data = np.array(a.get(y_gene_name))

        x_mrna_data = np.array(b.get(x_gene_name))
        y_mrna_data = np.array(b.get(y_gene_name))

        old_data_list = [x_protein_data,y_protein_data,x_mrna_data,y_mrna_data]
        new_data_list = []
        for data in old_data_list:
            if not np.all(data):
                data = np.zeros(122)
            new_data_list.append(data)

        krug_cds[4].data = {subtype_tuple: krug_subtype_tumor_tuple,
                                 subtype: list(zip(*krug_subtype_tumor_tuple))[0],
                                 'x_protein_data': new_data_list[0],
                                 'y_protein_data': new_data_list[1]}
        krug_cds[5].data = {subtype_tuple: krug_subtype_tumor_tuple,
                                 subtype: list(zip(*krug_subtype_tumor_tuple))[0],
                                 'x_mRNA_data': new_data_list[2],
                                 'y_mRNA_data': new_data_list[3]}

    with h5py.File(MertinsProteome, "r") as a, h5py.File(MertinsTranscriptome, "r") as b:
        x_protein_data = np.array(a.get(x_gene_name))
        y_protein_data = np.array(a.get(y_gene_name))

        x_mrna_data = np.array(b.get(x_gene_name))
        y_mrna_data = np.array(b.get(y_gene_name))

        old_data_list = [x_protein_data,y_protein_data,x_mrna_data,y_mrna_data]
        new_data_list = []
        for data in old_data_list:
            if not np.all(data):
                data = np.zeros(77)
            new_data_list.append(data)

        mertins_cds[4].data = {subtype_tuple: mertins_subtype_tumor_tuple,
                                 subtype: list(zip(*mertins_subtype_tumor_tuple))[0],
                                 'x_protein_data': new_data_list[0],
                                 'y_protein_data': new_data_list[1]}
        mertins_cds[5].data = {subtype_tuple: mertins_subtype_tumor_tuple,
                                 subtype: list(zip(*mertins_subtype_tumor_tuple))[0],
                                 'x_mRNA_data': new_data_list[2],
                                 'y_mRNA_data': new_data_list[3]}
    # scatter plot title change
    jo_pro_pro_scatter_plot.title.text = f"{y_gene_name} vs. {x_gene_name}"
    kr_pro_pro_scatter_plot.title.text = f"{y_gene_name} vs. {x_gene_name}"
    me_pro_pro_scatter_plot.title.text = f"{y_gene_name} vs. {x_gene_name}"

    jo_rna_rna_scatter_plot.title.text = f"{y_gene_name} vs. {x_gene_name}"
    kr_rna_rna_scatter_plot.title.text = f"{y_gene_name} vs. {x_gene_name}"
    me_rna_rna_scatter_plot.title.text = f"{y_gene_name} vs. {x_gene_name}"


def Johansson_Subtype_Plot_update():
    """ Create ColumnDataSource dictionary for plotting subtype mean and SEM purposes """
    # Create empty list to store data iterated from the for loop below
    x_axis_subtype = []

    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    for gene in SUBTYPE_PLOT_GENE_LIST:
        try:
            gene_sub_avg_pro = jo_pro.loc[gene,['basal_avg','her2_avg','lumA_avg','lumB_avg','norm_avg']].tolist()
            gene_sub_upper_pro = jo_pro.loc[gene, ['basal_up', 'her2_up', 'lumA_up', 'lumB_up', 'norm_up']].tolist()
            gene_sub_lower_pro = jo_pro.loc[gene, ['basal_down', 'her2_down', 'lumA_down', 'lumB_down', 'norm_down']].tolist()
        except KeyError:
            gene_sub_avg_pro = [0, 0, 0, 0, 0]
            gene_sub_upper_pro = [0, 0, 0, 0, 0]
            gene_sub_lower_pro = [0, 0, 0, 0, 0]

        try:
            gene_sub_avg_rna = jo_rna.loc[gene,['basal_avg','her2_avg','lumA_avg','lumB_avg','norm_avg']].tolist()
            gene_sub_upper_rna = jo_rna.loc[gene, ['basal_up', 'her2_up', 'lumA_up', 'lumB_up', 'norm_up']].tolist()
            gene_sub_lower_rna = jo_rna.loc[gene, ['basal_down', 'her2_down', 'lumA_down', 'lumB_down', 'norm_down']].tolist()
        except KeyError:
            gene_sub_avg_rna = [0, 0, 0, 0, 0]
            gene_sub_upper_rna = [0, 0, 0, 0, 0]
            gene_sub_lower_rna = [0, 0, 0, 0, 0]

        x_axis_subtype.extend([(gene,x) for x in ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']])

        y_axis_subtype_protein.extend(gene_sub_avg_pro)
        y_axis_subtype_mrna.extend(gene_sub_avg_rna)

        upper_bar_protein.extend(float(x) for x in gene_sub_upper_pro)
        lower_bar_protein.extend(float(x) for x in gene_sub_lower_pro)

        upper_bar_mrna.extend(float(x) for x in gene_sub_upper_rna)
        lower_bar_mrna.extend(float(x) for x in gene_sub_lower_rna)

    # modify the plot data sources from the dictionaries for each gene.
    jo_protein_subtype_plot.x_range.factors = x_axis_subtype
    jo_subtype_protein_CDS.data = {'x': x_axis_subtype, 'y': y_axis_subtype_protein,
              'upper': upper_bar_protein, 'lower': lower_bar_protein}

    jo_mRNA_subtype_plot.x_range.factors = x_axis_subtype
    jo_subtype_mRNA_CDS.data = {'x': x_axis_subtype, 'y': y_axis_subtype_mrna,
              'upper': upper_bar_mrna, 'lower': lower_bar_mrna}

def Krug_Subtype_Plot_update():
    """ Create ColumnDataSource dictionary for plotting subtype mean and SEM purposes """
    # Create empty list to store data iterated from the for loop below
    x_axis_subtype = []

    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    for gene in SUBTYPE_PLOT_GENE_LIST:
        try:
            gene_sub_avg_pro = kr_pro.loc[gene,['basal_avg','her2_avg','lumA_avg','lumB_avg','norm_avg']].tolist()
            gene_sub_upper_pro = kr_pro.loc[gene,['basal_up','her2_up','lumA_up','lumB_up','norm_up']].tolist()
            gene_sub_lower_pro = kr_pro.loc[gene,['basal_down','her2_down','lumA_down','lumB_down','norm_down']].tolist()
        except KeyError:
            gene_sub_avg_pro = [0, 0, 0, 0, 0]
            gene_sub_upper_pro = [0, 0, 0, 0, 0]
            gene_sub_lower_pro = [0, 0, 0, 0, 0]

        try:
            gene_sub_avg_rna = kr_rna.loc[gene,['basal_avg','her2_avg','lumA_avg','lumB_avg','norm_avg']].tolist()
            gene_sub_upper_rna = kr_rna.loc[gene,['basal_up','her2_up','lumA_up','lumB_up','norm_up']].tolist()
            gene_sub_lower_rna = kr_rna.loc[gene,['basal_down','her2_down','lumA_down','lumB_down','norm_down']].tolist()

        except KeyError:
            gene_sub_avg_rna = [0, 0, 0, 0, 0]
            gene_sub_upper_rna = [0, 0, 0, 0, 0]
            gene_sub_lower_rna = [0, 0, 0, 0, 0]

        x_axis_subtype.extend([(gene,x) for x in ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']])

        y_axis_subtype_protein.extend(gene_sub_avg_pro)
        y_axis_subtype_mrna.extend(gene_sub_avg_rna)

        upper_bar_protein.extend(float(x) for x in gene_sub_upper_pro)
        lower_bar_protein.extend(float(x) for x in gene_sub_lower_pro)

        upper_bar_mrna.extend(float(x) for x in gene_sub_upper_rna)
        lower_bar_mrna.extend(float(x) for x in gene_sub_lower_rna)

    # modify the plot data sources from the dictionaries for each gene.
    kr_protein_subtype_plot.x_range.factors = x_axis_subtype
    kr_subtype_protein_CDS.data = {'x': x_axis_subtype, 'y': y_axis_subtype_protein,
                                   'upper': upper_bar_protein, 'lower': lower_bar_protein}

    kr_mRNA_subtype_plot.x_range.factors = x_axis_subtype
    kr_subtype_mRNA_CDS.data = {'x': x_axis_subtype, 'y': y_axis_subtype_mrna,
                                'upper': upper_bar_mrna, 'lower': lower_bar_mrna}

def Mertins_Subtype_Plot_update():
    """ Create ColumnDataSource dictionary for plotting subtype mean and SEM purposes """
    # Create empty list to store data iterated from the for loop below
    x_axis_subtype = []

    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    for gene in SUBTYPE_PLOT_GENE_LIST:
        try:
            gene_sub_avg_pro = me_pro.loc[gene,['basal_avg','her2_avg','lumA_avg','lumB_avg']].tolist()
            gene_sub_upper_pro = me_pro.loc[gene,['basal_up','her2_up','lumA_up','lumB_up']].tolist()
            gene_sub_lower_pro = me_pro.loc[gene,['basal_down','her2_down','lumA_down','lumB_down']].tolist()

        except KeyError:
            gene_sub_avg_pro = [0, 0, 0, 0]
            gene_sub_upper_pro = [0, 0, 0, 0]
            gene_sub_lower_pro = [0, 0, 0, 0]

        try:
            gene_sub_avg_rna = me_rna.loc[gene,['basal_avg','her2_avg','lumA_avg','lumB_avg']].tolist()
            gene_sub_upper_rna = me_rna.loc[gene,['basal_up','her2_up','lumA_up','lumB_up']].tolist()
            gene_sub_lower_rna = me_rna.loc[gene,['basal_down','her2_down','lumA_down','lumB_down']].tolist()

        except KeyError:
            gene_sub_avg_rna = [0, 0, 0, 0]
            gene_sub_upper_rna = [0, 0, 0, 0]
            gene_sub_lower_rna = [0, 0, 0, 0]

        x_axis_subtype.extend([(gene,x) for x in ['Basal', 'Her2', 'LumA', 'LumB']])

        y_axis_subtype_protein.extend(gene_sub_avg_pro)
        y_axis_subtype_mrna.extend(gene_sub_avg_rna)

        upper_bar_protein.extend(float(x) for x in gene_sub_upper_pro)
        lower_bar_protein.extend(float(x) for x in gene_sub_lower_pro)

        upper_bar_mrna.extend(float(x) for x in gene_sub_upper_rna)
        lower_bar_mrna.extend(float(x) for x in gene_sub_lower_rna)

    # modify the plot data sources from the dictionaries for each gene.
    me_protein_subtype_plot.x_range.factors = x_axis_subtype
    me_subtype_protein_CDS.data = {'x': x_axis_subtype, 'y': y_axis_subtype_protein,
                                   'upper': upper_bar_protein, 'lower': lower_bar_protein}

    me_mRNA_subtype_plot.x_range.factors = x_axis_subtype
    me_subtype_mRNA_CDS.data = {'x': x_axis_subtype, 'y': y_axis_subtype_mrna,
                                'upper': upper_bar_mrna, 'lower': lower_bar_mrna}

def Subtype_Plot_GeneList_Update(attrname, old, new):
    global SUBTYPE_PLOT_GENE_LIST
    SUBTYPE_PLOT_GENE_LIST[TICKER_INDEX] = new
    print(f'Subtype gene list updated: {SUBTYPE_PLOT_GENE_LIST}...')

def Master_Subtype_Plot_Button_Update(event):
    start = time.time()
    Johansson_Subtype_Plot_update()
    Krug_Subtype_Plot_update()
    Mertins_Subtype_Plot_update()
    end = time.time()
    print(f"Subtype plot updated. {round(end - start, 6)} ms")

# --------------------------------------------------------------------
# Master update that for text boxes
# --------------------------------------------------------------------
def Master_Textbox_Update(attrname, old, new):
    start = time.time()

    Line_Plot_Update(attrname, old, new)
    #Scatter_Plot_Update(attrname, old, new)
    Subtype_Plot_GeneList_Update(attrname, old, new)
    Textbox_GeneList_Update(attrname, old, new)
    STATUS_FUNCTION_LIST[TICKER_INDEX](attrname, old, new)

    end = time.time()
    print(f"It takes {round(end-start, 6)} ms")

# --------------------------------------------------------------------
# Correlation tables
# --------------------------------------------------------------------
def Correlation_Table_Protein_Update(attrname, old, new):
    start = time.time()
    [jo_Correlation_table_df, kr_Correlation_table_df, me_Correlation_table_df, protein_sum_df] = Get_Protein_Correlation_Table(new)
    jo_pro_gene_cor_source.data = jo_Correlation_table_df # update column data source
    jo_pro_cor_data_table.source = jo_pro_gene_cor_source # update table

    kr_pro_gene_cor_source.data = kr_Correlation_table_df  # update column data source
    kr_pro_cor_data_table.source = kr_pro_gene_cor_source  # update table

    me_pro_gene_cor_source.data = me_Correlation_table_df  # update column data source
    me_pro_cor_data_table.source = me_pro_gene_cor_source  # update table

    summary_pro_source.data = protein_sum_df  # update column data source
    summary_pro_data_table.source = summary_pro_source  # update summary table

    #protein_table_gene_col_name.title = f'Genes correlated with {new}' # update plot title
    end = time.time()
    print(f"Correlation Table Protein is updated. {round(end-start, 6)} ms")

def Correlation_Table_mRNA_Update(attrname, old, new):
    start = time.time()
    [jo_m_Correlation_table_df, kr_m_Correlation_table_df, me_m_Correlation_table_df, mRNA_sum_df] = Get_mRNA_Correlation_Table(new)
    jo_mrna_gene_cor_source.data = jo_m_Correlation_table_df # update column data source
    jo_mrna_cor_data_table.source = jo_mrna_gene_cor_source # update table

    kr_mrna_gene_cor_source.data = kr_m_Correlation_table_df  # update column data source
    kr_mrna_cor_data_table.source = kr_mrna_gene_cor_source  # update table

    me_mrna_gene_cor_source.data = me_m_Correlation_table_df  # update column data source
    me_mrna_cor_data_table.source = me_mrna_gene_cor_source  # update table

    summary_mrna_source.data = mRNA_sum_df  # update column data source
    summary_mrna_data_table.source = summary_mrna_source  # update summary table

    #mRNA_table_gene_col_name.title = f'Genes correlated with {new}' # update plot title
    end = time.time()
    print(f"Correlation Table mRNA is updated. {round(end-start, 6)} ms")

# --------------------------------------------------------------------
# Row 1: Protein Complex Subunit Correlation Plot
# --------------------------------------------------------------------
start_time = time.time()

# function call
[jo_plot_p, jo_plot_m] = Johansson_Line_Plot()
[kr_plot_p, kr_plot_m] = Krug_Line_Plot()
[me_plot_p, me_plot_m] = Mertins_Line_Plot()

[johansson_cds, johansson_subtype_tumor_tuple, jo_unique_gene_list] = Johansson_CDS()
[krug_cds, krug_subtype_tumor_tuple, kr_unique_gene_list] = Krug_CDS()
[mertins_cds, mertins_subtype_tumor_tuple, me_unique_gene_list] = Mertins_CDS()

all_unique_genes = Gene_List()

# constants
GENE_NUMBER = ['Gene 1', 'Gene 2', 'Gene 3', 'Gene 4']
INITIAL_GENE = ['NDUFS2', 'NDUFS3', 'NDUFS7', 'NDUFS8']
GENE_COLORS = ['red', 'blue', 'green', 'orange']

for j in range(4):
    jo_plot_p.line(subtype_tuple, protein_data, source=johansson_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    jo_plot_m.line(subtype_tuple, mRNA_data, source=johansson_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    TICKER[j] = AutocompleteInput(title=GENE_NUMBER[j], value=INITIAL_GENE[j], width=100, width_policy='fixed',
                                  min_characters=2, completions=Nix(INITIAL_GENE[j], all_unique_genes),
                                  case_sensitive=False)
    TICKER[j].on_change('value', TICKER_FUNCTION_LIST[j])
    TICKER[j].on_change('value', Master_Textbox_Update)


for j in range(4):
    kr_plot_p.line(subtype_tuple, protein_data, source=krug_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    kr_plot_m.line(subtype_tuple, mRNA_data, source=krug_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    TICKER[j] = AutocompleteInput(title=GENE_NUMBER[j], value=INITIAL_GENE[j], width=100, width_policy='fixed',
                                  min_characters=2, completions=Nix(INITIAL_GENE[j], all_unique_genes),
                                  case_sensitive=False)
    TICKER[j].on_change('value', TICKER_FUNCTION_LIST[j])
    TICKER[j].on_change('value', Master_Textbox_Update)

for j in range(4):
    me_plot_p.line(subtype_tuple, protein_data, source=mertins_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    me_plot_m.line(subtype_tuple, mRNA_data, source=mertins_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    TICKER[j] = AutocompleteInput(title=GENE_NUMBER[j], value=INITIAL_GENE[j], width=100, width_policy='fixed',
                                  min_characters=2, completions=Nix(INITIAL_GENE[j], all_unique_genes),
                                  case_sensitive=False)
    TICKER[j].on_change('value', TICKER_FUNCTION_LIST[j])
    TICKER[j].on_change('value', Master_Textbox_Update)

# styling plots
for i in [jo_plot_p, jo_plot_m, kr_plot_p, kr_plot_m, me_plot_p, me_plot_m]:
    StylePlots(i, PlotID='Correlation')

# Put protein plot and mRNA plot into the same column, with a spacer separating them
line_plot_jo_layout = layout(column(jo_plot_p, Spacer(height=30), jo_plot_m))
line_plot_kr_layout = layout(column(kr_plot_p, Spacer(height=30), kr_plot_m))
line_plot_me_layout = layout(column(me_plot_p, Spacer(height=30), me_plot_m))

# Stack the 3 layouts into tabs
line_plot_tab = Tabs(tabs=[TabPanel(child=line_plot_jo_layout, title="Johansson"),
                     TabPanel(child=line_plot_kr_layout, title="Krug"),
                     TabPanel(child=line_plot_me_layout, title="Mertins")])

#Download data
button1 = Button(label="Download Data", button_type="success", width=120, width_policy='fixed')  # Gene 1 Data
button1.js_on_event("button_click", CustomJS(args=dict(source=johansson_cds[0]), code=open(join(dirname(__file__),
                                                                                                "Download_Javascript/jo_download.js")).read()))
button1.js_on_event("button_click", CustomJS(args=dict(source=krug_cds[0]), code=open(join(dirname(__file__),
                                                                                           "Download_Javascript/kr_download.js")).read()))
button1.js_on_event("button_click", CustomJS(args=dict(source=mertins_cds[0]), code=open(join(dirname(__file__),
                                                                                              "Download_Javascript/me_download.js")).read()))
####
button2 = Button(label="Download Data", button_type="success", width=120, width_policy='fixed')  # Gene 2 Data
button2.js_on_event("button_click", CustomJS(args=dict(source=johansson_cds[1]), code=open(join(dirname(__file__),
                                                                                                "Download_Javascript/jo_download.js")).read()))
button2.js_on_event("button_click", CustomJS(args=dict(source=krug_cds[1]), code=open(join(dirname(__file__),
                                                                                           "Download_Javascript/kr_download.js")).read()))
button2.js_on_event("button_click", CustomJS(args=dict(source=mertins_cds[1]), code=open(join(dirname(__file__),
                                                                                              "Download_Javascript/me_download.js")).read()))
####
button3 = Button(label="Download Data", button_type="success", width=120, width_policy='fixed')  # Gene 3 Data
button3.js_on_event("button_click", CustomJS(args=dict(source=johansson_cds[2]), code=open(join(dirname(__file__),
                                                                                                "Download_Javascript/jo_download.js")).read()))
button3.js_on_event("button_click", CustomJS(args=dict(source=krug_cds[2]), code=open(join(dirname(__file__),
                                                                                           "Download_Javascript/kr_download.js")).read()))
button3.js_on_event("button_click", CustomJS(args=dict(source=mertins_cds[2]), code=open(join(dirname(__file__),
                                                                                              "Download_Javascript/me_download.js")).read()))
####
button4 = Button(label="Download Data", button_type="success", width=120, width_policy='fixed')  # Gene 4 Data
button4.js_on_event("button_click", CustomJS(args=dict(source=johansson_cds[3]), code=open(join(dirname(__file__),
                                                                                                "Download_Javascript/jo_download.js")).read()))
button4.js_on_event("button_click", CustomJS(args=dict(source=krug_cds[3]), code=open(join(dirname(__file__),
                                                                                           "Download_Javascript/kr_download.js")).read()))
button4.js_on_event("button_click", CustomJS(args=dict(source=mertins_cds[3]), code=open(join(dirname(__file__),
                                                                                              "Download_Javascript/me_download.js")).read()))

# textbox and download buttons next to line plot --> TextBox_and_buttons
formatted_button1 = row(column(Spacer(height=19),button1))
formatted_button2 = row(column(Spacer(height=19),button2))
formatted_button3 = row(column(Spacer(height=19),button3))
formatted_button4 = row(column(Spacer(height=19),button4))

tickers_buttons_layout = column(row(TICKER[0], formatted_button1), row(TICKER[1], formatted_button2), row(TICKER[2], formatted_button3), row(TICKER[3], formatted_button4))

# --------------------------------------------------------------------
# Row 1.5: Pro-Pro & RNA-RNA Correlation Scatter Plot
# --------------------------------------------------------------------
# function call to get the plots
[jo_pro_pro_scatter_plot, kr_pro_pro_scatter_plot, me_pro_pro_scatter_plot] = Pro_Pro_Scatter_Plot()
[jo_rna_rna_scatter_plot, kr_rna_rna_scatter_plot, me_rna_rna_scatter_plot] = RNA_RNA_Scatter_Plot()

# constants
five_subtypes = ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']  # list of subtypes
five_subtype_colors = ['#E31A1C', '#FB9A99', '#1F78B4', '#A6CEE3','#33A02C']

# scatter plots
jo_pro_pro_scatter_plot.scatter('x_protein_data', 'y_protein_data', source=johansson_cds[4],
                                color=factor_cmap(subtype, five_subtype_colors, five_subtypes))
kr_pro_pro_scatter_plot.scatter('x_protein_data', 'y_protein_data', source=krug_cds[4],
                                color=factor_cmap(subtype, five_subtype_colors, five_subtypes))
me_pro_pro_scatter_plot.scatter('x_protein_data', 'y_protein_data', source=mertins_cds[4],
                                color=factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]))

jo_rna_rna_scatter_plot.scatter('x_mRNA_data', 'y_mRNA_data', source=johansson_cds[5],
                                color=factor_cmap(subtype, five_subtype_colors, five_subtypes))
kr_rna_rna_scatter_plot.scatter('x_mRNA_data', 'y_mRNA_data', source=krug_cds[5],
                                color=factor_cmap(subtype, five_subtype_colors, five_subtypes))
me_rna_rna_scatter_plot.scatter('x_mRNA_data', 'y_mRNA_data', source=mertins_cds[5],
                                color=factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]))

# put plots into a list and run through styling function
p_p_plots = [jo_pro_pro_scatter_plot, kr_pro_pro_scatter_plot, me_pro_pro_scatter_plot]
r_r_plots = [jo_rna_rna_scatter_plot, kr_rna_rna_scatter_plot, me_rna_rna_scatter_plot]
for i,j in zip(p_p_plots,r_r_plots):
    StylePlots(i, PlotID='mRNA-Prot')
    StylePlots(j, PlotID='mRNA-Prot')

# layout
scatter_plot_jo_layout = layout(column(row(p_p_plots[0]), Spacer(height=10),row(r_r_plots[0])))

scatter_plot_kr_layout = layout(column(row(p_p_plots[1]), Spacer(height=10),row(r_r_plots[1])))

scatter_plot_me_layout = layout(column(row(p_p_plots[2]), Spacer(height=10),row(r_r_plots[2])))

# stack 3 layouts into tab
pprr_scatter_plot_tab = Tabs(tabs=[TabPanel(child=scatter_plot_jo_layout, title="Johansson"),
                              TabPanel(child=scatter_plot_kr_layout, title="Krug"),
                              TabPanel(child=scatter_plot_me_layout, title="Mertins")])

# select widget
x_axis_select = Select(title="X axis", value="Gene 1", options=['Gene 1', 'Gene 2', 'Gene 3', 'Gene 4'], width = 100, width_policy='fixed')
y_axis_select = Select(title="Y axis", value="Gene 2", options=['Gene 1', 'Gene 2', 'Gene 3', 'Gene 4'], width = 100, width_policy='fixed')
selection_update_button = Button(label='Update', button_type='danger', width=100, width_policy='fixed')
selection_update_button.on_click(Scatter_Select_Button_Update)

# select widget next to pro-pro/rna-rna scatter plot
select_widget_layout = column(x_axis_select, y_axis_select, selection_update_button)

# static legend under selection dropdown
legend_pic = Div(text="<img src='https://i.imgur.com/bHTfVP8.png'>")

# --------------------------------------------------------------------
# Row 2: mRNA-Protein Correlation Scatter Plot
# --------------------------------------------------------------------
# function calls
[jo_plot_mRNA_prot1, jo_plot_mRNA_prot2, jo_plot_mRNA_prot3, jo_plot_mRNA_prot4] = Johansson_RNA_Pro_Scatter_Plot()
[kr_plot_mRNA_prot1, kr_plot_mRNA_prot2, kr_plot_mRNA_prot3, kr_plot_mRNA_prot4] = Krug_RNA_Pro_Scatter_Plot()
[me_plot_mRNA_prot1, me_plot_mRNA_prot2, me_plot_mRNA_prot3, me_plot_mRNA_prot4] = Mertins_RNA_Pro_Scatter_Plot()

johansson_scatter_cds = Johansson_Scatter_CDS()
krug_scatter_cds = Krug_Scatter_CDS()
mertins_scatter_cds = Mertins_Scatter_CDS()

jo_plot_mRNA_prot = [jo_plot_mRNA_prot1, jo_plot_mRNA_prot2, jo_plot_mRNA_prot3, jo_plot_mRNA_prot4]
kr_plot_mRNA_prot = [kr_plot_mRNA_prot1, kr_plot_mRNA_prot2, kr_plot_mRNA_prot3, kr_plot_mRNA_prot4]
me_plot_mRNA_prot = [me_plot_mRNA_prot1, me_plot_mRNA_prot2, me_plot_mRNA_prot3, me_plot_mRNA_prot4]

# scatter plots
jo_plot_mRNA_prot1.scatter(mRNA_data, protein_data, source=johansson_scatter_cds[0],
                           color=factor_cmap(subtype, five_subtype_colors, five_subtypes))
jo_plot_mRNA_prot2.scatter(mRNA_data, protein_data, source=johansson_scatter_cds[1],
                           color=factor_cmap(subtype, five_subtype_colors, five_subtypes))
jo_plot_mRNA_prot3.scatter(mRNA_data, protein_data, source=johansson_scatter_cds[2],
                           color=factor_cmap(subtype, five_subtype_colors, five_subtypes))
jo_plot_mRNA_prot4.scatter(mRNA_data, protein_data, source=johansson_scatter_cds[3],
                           color=factor_cmap(subtype, five_subtype_colors, five_subtypes))


kr_plot_mRNA_prot1.scatter(mRNA_data, protein_data, source=krug_scatter_cds[0],
                           color=factor_cmap(subtype, five_subtype_colors, five_subtypes))
kr_plot_mRNA_prot2.scatter(mRNA_data, protein_data, source=krug_scatter_cds[1],
                           color=factor_cmap(subtype, five_subtype_colors, five_subtypes))
kr_plot_mRNA_prot3.scatter(mRNA_data, protein_data, source=krug_scatter_cds[2],
                           color=factor_cmap(subtype, five_subtype_colors, five_subtypes))
kr_plot_mRNA_prot4.scatter(mRNA_data, protein_data, source=krug_scatter_cds[3],
                           color=factor_cmap(subtype, five_subtype_colors, five_subtypes))


me_plot_mRNA_prot1.scatter(mRNA_data, protein_data, source=mertins_scatter_cds[0],
                           color=factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]))
me_plot_mRNA_prot2.scatter(mRNA_data, protein_data, source=mertins_scatter_cds[1],
                           color=factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]))
me_plot_mRNA_prot3.scatter(mRNA_data, protein_data, source=mertins_scatter_cds[2],
                           color=factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]))
me_plot_mRNA_prot4.scatter(mRNA_data, protein_data, source=mertins_scatter_cds[3],
                           color=factor_cmap(subtype, five_subtype_colors[0:4], five_subtypes[0:4]))

# styling plots
for i in [jo_plot_mRNA_prot1, jo_plot_mRNA_prot2, jo_plot_mRNA_prot3, jo_plot_mRNA_prot4]:
    StylePlots(i, PlotID='mRNA-Prot')
for i in [kr_plot_mRNA_prot1, kr_plot_mRNA_prot2, kr_plot_mRNA_prot3, kr_plot_mRNA_prot4]:
    StylePlots(i, PlotID='mRNA-Prot')
for i in [me_plot_mRNA_prot1, me_plot_mRNA_prot2, me_plot_mRNA_prot3, me_plot_mRNA_prot4]:
    StylePlots(i, PlotID='mRNA-Prot')

# put gene1 and gene2 mRNA-protein correlation plots in the same row, with a spacer between
# put gene3 and gene4 mRNA-protein correlation plots in the row below, with a spacer between.
# spacer to divide the rows
scatter_plot_jo_layout = layout(column(row(jo_plot_mRNA_prot1,jo_plot_mRNA_prot2),
                                       Spacer(height=10),
                                       row(jo_plot_mRNA_prot3,jo_plot_mRNA_prot4)))
scatter_plot_kr_layout = layout(column(row(kr_plot_mRNA_prot1,kr_plot_mRNA_prot2),
                                       Spacer(height=10),
                                       row(kr_plot_mRNA_prot3,kr_plot_mRNA_prot4)))
scatter_plot_me_layout = layout(column(row(me_plot_mRNA_prot1,me_plot_mRNA_prot2),
                                       Spacer(height=10),
                                       row(me_plot_mRNA_prot3,me_plot_mRNA_prot4)))
# stack 3 layouts into tab
scatter_plot_tab = Tabs(tabs=[TabPanel(child=scatter_plot_jo_layout, title="Johansson"),
                              TabPanel(child=scatter_plot_kr_layout, title="Krug"),
                              TabPanel(child=scatter_plot_me_layout, title="Mertins")])

# scatter plots update button
scatter_update_button = Button(label='Update', button_type='danger', width=50, width_policy='fixed')
scatter_update_button.on_click(Scatter_Plot_Update)

# --------------------------------------------------------------------
# Row 3: Subtype Mean and SEM plot
# --------------------------------------------------------------------
# function call
[jo_subtype_protein_CDS, jo_subtype_mRNA_CDS] = Johansson_Subtype_Avg_Plot_CDS()
[kr_subtype_protein_CDS, kr_subtype_mRNA_CDS] = Krug_Subtype_Avg_Plot_CDS()
[me_subtype_protein_CDS, me_subtype_mRNA_CDS] = Mertins_Subtype_Avg_Plot_CDS()

[jo_protein_subtype_plot, jo_mRNA_subtype_plot] = Johansson_Subtype_Plot(jo_subtype_protein_CDS, jo_subtype_mRNA_CDS)
[kr_protein_subtype_plot, kr_mRNA_subtype_plot] = Krug_Subtype_Plot(kr_subtype_protein_CDS, kr_subtype_mRNA_CDS)
[me_protein_subtype_plot, me_mRNA_subtype_plot] = Mertins_Subtype_Plot(me_subtype_protein_CDS, me_subtype_mRNA_CDS)

# set aesthetic values
jo_plot_subtype_p = StylePlots(jo_protein_subtype_plot, PlotID='Subtype')
jo_plot_subtype_m = StylePlots(jo_mRNA_subtype_plot, PlotID='Subtype')

kr_plot_subtype_p = StylePlots(kr_protein_subtype_plot, PlotID='Subtype')
kr_plot_subtype_m = StylePlots(kr_mRNA_subtype_plot, PlotID='Subtype')

me_plot_subtype_p = StylePlots(me_protein_subtype_plot, PlotID='Subtype')
me_plot_subtype_m = StylePlots(me_mRNA_subtype_plot, PlotID='Subtype')

# Put the subtype plots for protein and mRNA into the same row, each plot takes a column, with a spacer separating them
subtype_plot_jo_layout = layout(row(column(jo_plot_subtype_p), column(jo_plot_subtype_m)))
subtype_plot_kr_layout = layout(row(column(kr_plot_subtype_p), column(kr_plot_subtype_m)))
subtype_plot_me_layout = layout(row(column(me_plot_subtype_p), column(me_plot_subtype_m)))

# Stack the 3 layouts into tabs
subtype_plot_tab = Tabs(tabs=[TabPanel(child=subtype_plot_jo_layout, title="Johansson"),
                              TabPanel(child=subtype_plot_kr_layout, title="Krug"),
                              TabPanel(child=subtype_plot_me_layout, title="Mertins")])

# button to update subtype plot (beta)
subtype_update_button = Button(label='Update', button_type='danger', width=50, width_policy='fixed')
subtype_update_button.on_click(Master_Subtype_Plot_Button_Update)


# --------------------------------------------------------------------
# Row 4: GxG Correlation
# --------------------------------------------------------------------
# function call
[jo_mrna_ERBB2, jo_protein_ERBB2, kr_mrna_ERBB2, kr_protein_ERBB2, me_mrna_ERBB2, me_protein_ERBB2, summary_protein_ERBB2, summary_mrna_ERBB2] = Import_Static_Correlation_Table()

# constant
TABLE_WIDTH = 900
TABLE_HEIGHT = 250

# textbox widget and on_change
protein_correlation_textbox = AutocompleteInput(title='Gene:', value=str('ERBB2'), width=180, width_policy='fixed',
                                            min_characters=2, completions = all_unique_genes, case_sensitive=False)
protein_correlation_textbox.on_change('value', Correlation_Table_Protein_Update)
mrna_correlation_textbox = AutocompleteInput(title='Gene:', value=str('ERBB2'), width=180, width_policy='fixed',
                                            min_characters=2, completions = all_unique_genes, case_sensitive=False)
mrna_correlation_textbox.on_change('value', Correlation_Table_mRNA_Update)

# CDS
jo_pro_gene_cor_source = ColumnDataSource(data=jo_protein_ERBB2) # set up columndatasource
kr_pro_gene_cor_source = ColumnDataSource(data=kr_protein_ERBB2)
me_pro_gene_cor_source = ColumnDataSource(data=me_protein_ERBB2)
summary_pro_source = ColumnDataSource(data=summary_protein_ERBB2)

jo_mrna_gene_cor_source = ColumnDataSource(data=jo_mrna_ERBB2) # set up columndatasource
kr_mrna_gene_cor_source = ColumnDataSource(data=kr_mrna_ERBB2)
me_mrna_gene_cor_source = ColumnDataSource(data=me_mrna_ERBB2)
summary_mrna_source = ColumnDataSource(data=summary_mrna_ERBB2)

# table assembly
# title = f'Genes correlated with {mrna_correlation_textbox.value}
protein_table_gene_col_name = TableColumn(field='Gene', title='Gene')
mRNA_table_gene_col_name = TableColumn(field='Gene', title='Gene')

protein_table_columns = [protein_table_gene_col_name,  # set up protein_table_columns
                         TableColumn(field='r', title='Coefficient (Pearson)', formatter=ScientificFormatter(precision=3)),
                         TableColumn(field='p', title='p value', formatter=ScientificFormatter(precision=3)),
                         TableColumn(field='Name', title='Name')]
mRNA_table_columns = [mRNA_table_gene_col_name,  # set up mRNA_table_columns
                      TableColumn(field='r', title='Coefficient (Pearson)', formatter=ScientificFormatter(precision=3)),
                      TableColumn(field='p', title='p value', formatter=ScientificFormatter(precision=3)),
                      TableColumn(field='Name', title='Name')]

protein_summary_table_columns = [protein_table_gene_col_name, # set up protein_summary_table_columns
                         TableColumn(field='p avg', title='Mean of 2 smallest p values', formatter=ScientificFormatter(precision=3)),
                         TableColumn(field='Name', title='Name')]

mRNA_summary_table_columns = [mRNA_table_gene_col_name, # set up mRNA_summary_table_columns
                         TableColumn(field='p avg', title='Mean of 2 smallest p values', formatter=ScientificFormatter(precision=3)),
                         TableColumn(field='Name', title='Name')]

# create table widget by assembling columndatasource and protein_table_columns
jo_pro_cor_data_table = DataTable(source=jo_pro_gene_cor_source, columns=protein_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, index_position=None, autosize_mode='fit_columns')
kr_pro_cor_data_table = DataTable(source=kr_pro_gene_cor_source, columns=protein_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, index_position=None, autosize_mode='fit_columns')
me_pro_cor_data_table = DataTable(source=me_pro_gene_cor_source, columns=protein_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, index_position=None, autosize_mode='fit_columns')
summary_pro_data_table = DataTable(source=summary_pro_source, columns=protein_summary_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, index_position=None, autosize_mode='fit_columns')

# create table widget by assembling columndatasource and mRNA_table_columns
jo_mrna_cor_data_table = DataTable(source=jo_mrna_gene_cor_source, columns=mRNA_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, index_position=None, autosize_mode='fit_columns')
kr_mrna_cor_data_table = DataTable(source=kr_mrna_gene_cor_source, columns=mRNA_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, index_position=None, autosize_mode='fit_columns')
me_mrna_cor_data_table = DataTable(source=me_mrna_gene_cor_source, columns=mRNA_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, index_position=None, autosize_mode='fit_columns')
summary_mrna_data_table = DataTable(source=summary_mrna_source, columns=mRNA_summary_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, index_position=None, autosize_mode='fit_columns')

# combine tables into tabs
cor_pro_table_tab = Tabs(tabs=[
                           TabPanel(child=summary_pro_data_table, title = 'Summary'),
                           TabPanel(child=jo_pro_cor_data_table, title ="Johansson"),
                           TabPanel(child=kr_pro_cor_data_table, title = "Krug"),
                           TabPanel(child=me_pro_cor_data_table, title = "Mertins")
                           ])

cor_mrna_table_tab = Tabs(tabs=[
                           TabPanel(child=summary_mrna_data_table, title = 'Summary'),
                           TabPanel(child=jo_mrna_cor_data_table, title ="Johansson"),
                           TabPanel(child=kr_mrna_cor_data_table, title = "Krug"),
                           TabPanel(child=me_mrna_cor_data_table, title = "Mertins")
                           ])

# button and call back - download button will be removed
button5 = Button(label="Download Data", button_type="success", width=100, width_policy='fixed') #Protein Correlation Table
button5.js_on_event("button_click", CustomJS(args=dict(source=summary_pro_source),
                                             code=open(join(dirname(__file__),
                                                            "Download_Javascript/pro_sum_cor_table_download.js")).read()))


button6 = Button(label="Download Data", button_type="success", width=100, width_policy='fixed') #mRNA Correlation Table
button6.js_on_event("button_click", CustomJS(args=dict(source=summary_mrna_source),
                                             code=open(join(dirname(__file__),
                                                            "Download_Javascript/mRNA_sum_cor_table_download.js")).read()))

formatted_button5 = row(column(Spacer(height=19),button5))
formatted_button6 = row(column(Spacer(height=19),button6))

# --------------------------------------------------------------------
# Row 5: Heatmaps (future feature)
# --------------------------------------------------------------------

# [jo_pro_hm, kr_pro_hm, me_pro_hm] = Get_Protein_Heatmaps()
# [jo_mRNA_hm, kr_mRNA_hm, me_mRNA_hm] = Get_mRNA_Heatmaps()
#
# protein_heatmap_tab = Tabs(tabs=[TabPanel(child=jo_pro_hm, title="Johansson"),
#                                TabPanel(child=kr_pro_hm, title="Krug"),
#                                TabPanel(child=me_pro_hm, title="Mertins")])
# mRNA_heatmap_tab = Tabs(tabs=[TabPanel(child=jo_mRNA_hm, title="Johansson"),
#                                TabPanel(child=kr_mRNA_hm, title="Krug"),
#                                TabPanel(child=me_mRNA_hm, title="Mertins")])

# --------------------------------------------------------------------
# layout
# --------------------------------------------------------------------
RowSpacer = Spacer(height=30)
PageMargin = Spacer(width=30)

tool_title_div = Div(text="Breast Cancer Proteome Portal",
                     styles={'font-size': '200%', 'color': 'black', 'font-style': 'normal', 'font-weight': 'bold'},
                     width=1000)
tab_click_reminder_text_div = Div(text="Toggle between studies using the tabs above each plot.",
                         styles={'font-size': '100%', 'color': 'black'}, width=1000)

line_plot_div = column(Div(text="Abundance Traces:",
                           styles={'font-size': '120%', 'color': 'black', 'font-style': 'italic',
                                       'font-weight': 'bold'}))

scatter_plot_div = column(Div(text="Protein-mRNA:",
                              styles={'font-size': '120%', 'color': 'black', 'font-style': 'italic',
                                      'font-weight': 'bold'}))

scatter_plot_div2 = column(Div(text="Protein-Protein/mRNA-mRNA:",
                              styles={'font-size': '120%', 'color': 'black', 'font-style': 'italic',
                                      'font-weight': 'bold'}))

subtype_plot_div = column(Div(text="Abundances by Subtype:",
                              styles={'font-size': '120%', 'color': 'black', 'font-style': 'italic', 'font-weight': 'bold'}))

protein_cor_table_title_div = column(Div(text="Correlation Table - Protein:",
                                    styles={'font-size': '120%', 'color': 'black', 'font-style': 'italic', 'font-weight': 'bold'}))
mRNA_cor_table_title_div = column(Div(text="Correlation Table - mRNA:",
                                 styles={'font-size': '120%', 'color': 'black', 'font-style': 'italic', 'font-weight': 'bold'}))
gene_hide_text_div = Div(text="To hide traces, click on interactive legend entries.",
                         styles={'font-size': '100%', 'color': 'red'}, width=200, width_policy='auto', height_policy='auto')

tool_title_description_section = column(tool_title_div, tab_click_reminder_text_div)
gene_entry_section = column(tickers_buttons_layout)
line_plot_section = column(line_plot_div, row(line_plot_tab, column(gene_entry_section, gene_hide_text_div)))
select_widget_section = column(select_widget_layout)
scatter_plot_section = column(row(scatter_plot_div, scatter_update_button, Spacer(width=400), scatter_plot_div2), row(scatter_plot_tab, Spacer(width=120), pprr_scatter_plot_tab, column(select_widget_section, legend_pic)))
subtype_plot_section = row(column(row(subtype_plot_div, subtype_update_button), subtype_plot_tab),column(Spacer(height=60),legend_pic))
correlation_table_section = column(row(column(protein_cor_table_title_div, row(protein_correlation_textbox, formatted_button5), cor_pro_table_tab, Spacer(height=30),
                                              mRNA_cor_table_title_div, row(mrna_correlation_textbox, formatted_button6), cor_mrna_table_tab, Spacer(height=15))))

l = layout([
    [PageMargin, tool_title_description_section],
    Spacer(height=15),
    [PageMargin, line_plot_section],
    [RowSpacer],
    [PageMargin, scatter_plot_section],
    [RowSpacer],
    [PageMargin, subtype_plot_section],
    [RowSpacer],
    [PageMargin, correlation_table_section],
    [RowSpacer]])

curdoc().add_root(l)
curdoc().title="Breast Cancer Data Portal"

end_time=time.time()
elapsed=end_time-start_time
print(f'Total response time: {elapsed}')

# --------------------------------------------------------------------
# HTML related
# --------------------------------------------------------------------
