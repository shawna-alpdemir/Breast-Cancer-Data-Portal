##### Import ###########################################################################################################
import h5py
import numpy as np
import time
from os.path import join, dirname
from bokeh.io import curdoc
from bokeh.layouts import layout, column, row
from bokeh.models import Spacer, Tabs, Panel, AutocompleteInput, ColumnDataSource, TableColumn, DataTable, Button, \
    Select, Div, CustomJS, ScientificFormatter, PreText
from bokeh.transform import factor_cmap

from Correlation_Table import Get_Protein_Correlation_Table, Get_mRNA_Correlation_Table
from Gene_List import Gene_List
from Import_Files import Import_HDF5, Import_Static_Correlation_Table
from Plot_All_Styling import StylePlots
from Plot_Line_Scatter_ColumnDataSource import Johansson_CDS, Krug_CDS, Mertins_CDS, Johansson_Scatter_CDS, Krug_Scatter_CDS, Mertins_Scatter_CDS
from Line_Plot import Johansson_Line_Plot, Krug_Line_Plot, Mertins_Line_Plot
from Scatter_Plot import Pro_Pro_Scatter_Plot, RNA_RNA_Scatter_Plot, Johansson_RNA_Pro_Scatter_Plot, Krug_RNA_Pro_Scatter_Plot, Mertins_RNA_Pro_Scatter_Plot
from Subtype_Average_DF import Johansson_Subtype_Avg_SEM_DFs, Krug_Subtype_Avg_SEM_DFs, Mertins_Subtype_Avg_SEM_DFs
from Subtype_Average_Plot_CDS import Johansson_Subtype_Avg_Plot_CDS, Krug_Subtype_Avg_Plot_CDS, Mertins_Subtype_Avg_Plot_CDS
from Subtype_Plot import Johansson_Subtype_Plot, Krug_Subtype_Plot, Mertins_Subtype_Plot

##### Row 0: Functions Def ##############################3##############################################################
# constants
TICKER = {}
TICKER_INDEX = None
TICKER_GENE_LIST = ['ESR1', 'PGR', 'ERBB2', 'MKI67']

# dictionary key name
subtype_tuple = 'subtype, tumor'
subtype = 'subtype'
protein_data = 'protein_data'
mRNA_data = 'mRNA_data'
gene = 'gene'

# function call
[JohanssonProteome, JohanssonTranscriptome, KrugProteome, KrugTranscriptome, MertinsProteome, MertinsTranscriptome] = Import_HDF5()

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
    gene_hide_text_div.text = f"Gene 1 has been changed to {new}"

def Ticker2_Change(attrname, old, new):
    global TICKER_INDEX
    TICKER_INDEX = 1
    TICKER[0].completions = Nix(new, all_unique_genes)
    TICKER[2].completions = Nix(new, all_unique_genes)
    TICKER[3].completions = Nix(new, all_unique_genes)
    gene_hide_text_div.text = f"Gene 2 has been changed to {new}"

def Ticker3_Change(attrname, old, new):
    global TICKER_INDEX
    TICKER_INDEX = 2
    TICKER[0].completions = Nix(new, all_unique_genes)
    TICKER[1].completions = Nix(new, all_unique_genes)
    TICKER[3].completions = Nix(new, all_unique_genes)
    gene_hide_text_div.text = f"Gene 3 has been changed to {new}"

def Ticker4_Change(attrname, old, new):
    global TICKER_INDEX
    TICKER_INDEX = 3
    TICKER[0].completions = Nix(new, all_unique_genes)
    TICKER[1].completions = Nix(new, all_unique_genes)
    TICKER[2].completions = Nix(new, all_unique_genes)
    gene_hide_text_div.text = f"Gene 4 has been changed to {new}"

TICKER_FUNCTION_LIST = [Ticker1_Change, Ticker2_Change, Ticker3_Change, Ticker4_Change]

def Line_Plot_Update(attrname, old, new):
    """Update function for text entries call back targeting line plot"""
    start = time.time()

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

    end = time.time()
    print(f"line plots updated: {end-start}")

def Scatter_Plot_Update(attrname, old, new):
    """Update function for text entries call back targeting line plot"""
    start = time.time()

    with h5py.File(JohanssonProteome, "r") as a, h5py.File(JohanssonTranscriptome, "r") as b:
        protein_data_new_gene = np.array(a.get(new))
        mRNA_data_new_gene = np.array(b.get(new))

        if not np.all(protein_data_new_gene):
            protein_data_new_gene = np.zeros(45)

        if not np.all(mRNA_data_new_gene):
            mRNA_data_new_gene = np.zeroes(45)

        johansson_scatter_cds[TICKER_INDEX].data = {subtype: list(zip(*johansson_subtype_tumor_tuple))[0], # extract first element of the list of subtype_tumor_tuple, which is the subtype
                 protein_data: protein_data_new_gene,
                 mRNA_data: mRNA_data_new_gene}

    with h5py.File(KrugProteome, "r") as a, h5py.File(KrugTranscriptome, "r") as b:
        protein_data_new_gene = np.array(a.get(new))
        mRNA_data_new_gene = np.array(b.get(new))

        if not np.all(protein_data_new_gene):
            protein_data_new_gene = np.zeros(122)

        if not np.all(mRNA_data_new_gene):
            mRNA_data_new_gene = np.zeros(122)

        krug_scatter_cds[TICKER_INDEX].data = {subtype: list(zip(*krug_subtype_tumor_tuple))[0], # extract first element of the list of subtype_tumor_tuple
                 protein_data: protein_data_new_gene,
                 mRNA_data: mRNA_data_new_gene}

    with h5py.File(MertinsProteome, "r") as a, h5py.File(MertinsTranscriptome, "r") as b:
        protein_data_new_gene = np.array(a.get(new))
        mRNA_data_new_gene = np.array(b.get(new))

        if not np.all(protein_data_new_gene):
            protein_data_new_gene = np.zeros(77)

        if not np.all(mRNA_data_new_gene):
            mRNA_data_new_gene = np.zeros(77)

        mertins_scatter_cds[TICKER_INDEX].data = {subtype: list(zip(*mertins_subtype_tumor_tuple))[0], # extract first element of the list of subtype_tumor_tuple
                 protein_data: protein_data_new_gene,
                 mRNA_data: mRNA_data_new_gene}

    jo_plot_mRNA_prot[TICKER_INDEX].title.text = new
    kr_plot_mRNA_prot[TICKER_INDEX].title.text = new
    me_plot_mRNA_prot[TICKER_INDEX].title.text = new

    end = time.time()
    print(f"scatter plot updated {end-start}")

def Scatter_Select_Update(event):
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
                #data = np.repeat(np.nan, 45)
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
                #data = np.repeat(np.nan, 122)
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
                #data = np.repeat(np.nan, 77)
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
    jo_pro_pro_scatter_plot.title.text = f"{x_gene_name}-{y_gene_name}"
    kr_pro_pro_scatter_plot.title.text = f"{x_gene_name}-{y_gene_name}"
    me_pro_pro_scatter_plot.title.text = f"{x_gene_name}-{y_gene_name}"

    jo_rna_rna_scatter_plot.title.text = f"{x_gene_name}-{y_gene_name}"
    kr_rna_rna_scatter_plot.title.text = f"{x_gene_name}-{y_gene_name}"
    me_rna_rna_scatter_plot.title.text = f"{x_gene_name}-{y_gene_name}"

def Johansson_Subtype_Plot_update():
    x_axis_subtype = []

    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    jo_subtype = []

    for gene in TICKER_GENE_LIST:
        [jo_avg_protein_DF, jo_sem_protein_DF, jo_avg_mRNA_DF, jo_sem_mRNA_DF] = Johansson_Subtype_Avg_SEM_DFs(gene)

        x_axis_subtype.extend([(gene,x) for x in ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']])

        y_axis_subtype_protein.extend(jo_avg_protein_DF['values'].tolist())
        y_axis_subtype_mrna.extend(jo_avg_mRNA_DF['values'].tolist())

        avg_plus_sem_protein = jo_avg_protein_DF['values'] + jo_sem_protein_DF['values']
        avg_minus_sem_protein = jo_avg_protein_DF['values'] - jo_sem_protein_DF['values']

        avg_plus_sem_mRNA = jo_avg_mRNA_DF['values'] + jo_sem_mRNA_DF['values']
        avg_minus_sem_mRNA = jo_avg_mRNA_DF['values'] - jo_sem_mRNA_DF['values']

        upper_bar_protein.extend(avg_plus_sem_protein.tolist())
        lower_bar_protein.extend(avg_minus_sem_protein.to_list())

        upper_bar_mrna.extend(avg_plus_sem_mRNA.tolist())
        lower_bar_mrna.extend(avg_minus_sem_mRNA.to_list())

        jo_subtype.extend(['Basal', 'Her2', 'LumA', 'LumB', 'Norm'])

    # update plot and CDS
    jo_protein_subtype_plot.x_range.factors = x_axis_subtype
    jo_subtype_protein_CDS.data={'x': x_axis_subtype, 'y': y_axis_subtype_protein,
              'upper': upper_bar_protein, 'lower': lower_bar_protein}

    jo_mRNA_subtype_plot.x_range.factors = x_axis_subtype
    jo_subtype_mRNA_CDS.data = {'x': x_axis_subtype, 'y': y_axis_subtype_mrna,
                                       'upper': upper_bar_mrna, 'lower': lower_bar_mrna}

def Krug_Subtype_Plot_update():
    # Create empty list to store data iterated from the for loop below
    x_axis_subtype = []

    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    kr_subtype = []

    for gene in TICKER_GENE_LIST:
        [kr_avg_protein_DF, kr_sem_protein_DF, kr_avg_mRNA_DF, kr_sem_mRNA_DF] = Krug_Subtype_Avg_SEM_DFs(gene)

        x_axis_subtype.extend([(gene,x) for x in ['Basal', 'Her2', 'LumA', 'LumB', 'Norm']])

        y_axis_subtype_protein.extend(kr_avg_protein_DF['values'].tolist())
        y_axis_subtype_mrna.extend(kr_avg_mRNA_DF['values'].tolist())

        avg_plus_sem_protein = kr_avg_protein_DF['values'] + kr_sem_protein_DF['values']
        avg_minus_sem_protein = kr_avg_protein_DF['values'] - kr_sem_protein_DF['values']

        avg_plus_sem_mRNA = kr_avg_mRNA_DF['values'] + kr_sem_mRNA_DF['values']
        avg_minus_sem_mRNA = kr_avg_mRNA_DF['values'] - kr_sem_mRNA_DF['values']

        upper_bar_protein.extend(avg_plus_sem_protein.tolist())
        lower_bar_protein.extend(avg_minus_sem_protein.to_list())

        upper_bar_mrna.extend(avg_plus_sem_mRNA.tolist())
        lower_bar_mrna.extend(avg_minus_sem_mRNA.to_list())

        kr_subtype.extend(['Basal', 'Her2', 'LumA', 'LumB', 'Norm'])
    # update plot and CDS
    kr_protein_subtype_plot.x_range.factors = x_axis_subtype
    kr_subtype_protein_CDS.data={'x': x_axis_subtype, 'y': y_axis_subtype_protein,
              'upper': upper_bar_protein, 'lower': lower_bar_protein}

    kr_mRNA_subtype_plot.x_range.factors = x_axis_subtype
    kr_subtype_mRNA_CDS.data = {'x': x_axis_subtype, 'y': y_axis_subtype_mrna,
                                       'upper': upper_bar_mrna, 'lower': lower_bar_mrna}

def Mertins_Subtype_Plot_update():
    # Create empty list to store data iterated from the for loop below
    x_axis_subtype = []

    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    me_subtype = []

    for gene in TICKER_GENE_LIST:
        [me_avg_protein_DF, me_sem_protein_DF, me_avg_mRNA_DF, me_sem_mRNA_DF] = Mertins_Subtype_Avg_SEM_DFs(gene)

        x_axis_subtype.extend([(gene,x) for x in ['Basal', 'Her2', 'LumA', 'LumB']])

        y_axis_subtype_protein.extend(me_avg_protein_DF['values'].tolist())
        y_axis_subtype_mrna.extend(me_avg_mRNA_DF['values'].tolist())

        avg_plus_sem_protein = me_avg_protein_DF['values'] + me_sem_protein_DF['values']
        avg_minus_sem_protein = me_avg_protein_DF['values'] - me_sem_protein_DF['values']

        avg_plus_sem_mRNA = me_avg_mRNA_DF['values'] + me_sem_mRNA_DF['values']
        avg_minus_sem_mRNA = me_avg_mRNA_DF['values'] - me_sem_mRNA_DF['values']

        upper_bar_protein.extend(avg_plus_sem_protein.tolist())
        lower_bar_protein.extend(avg_minus_sem_protein.to_list())

        upper_bar_mrna.extend(avg_plus_sem_mRNA.tolist())
        lower_bar_mrna.extend(avg_minus_sem_mRNA.to_list())

        me_subtype.extend(['Basal', 'Her2', 'LumA', 'LumB'])
    # update plot and CDS
    me_protein_subtype_plot.x_range.factors = x_axis_subtype
    me_subtype_protein_CDS.data = {'x': x_axis_subtype, 'y': y_axis_subtype_protein,
                                   'upper': upper_bar_protein, 'lower': lower_bar_protein}

    me_mRNA_subtype_plot.x_range.factors = x_axis_subtype
    me_subtype_mRNA_CDS.data = {'x': x_axis_subtype, 'y': y_axis_subtype_mrna,
                                'upper': upper_bar_mrna, 'lower': lower_bar_mrna}

def All_Subtype_Plot_Update(attrname, old, new):
    start = time.time()
    global TICKER_GENE_LIST
    TICKER_GENE_LIST[TICKER_INDEX] = new
    Johansson_Subtype_Plot_update()
    Krug_Subtype_Plot_update()
    Mertins_Subtype_Plot_update()
    end = time.time()
    print(f'Current subtype plot gene list is: {TICKER_GENE_LIST}, {end-start}')

def Correlation_Table_Protein_Update(attrname, old, new):
    start = time.time()
    [jo_Correlation_table_df, kr_Correlation_table_df, me_Correlation_table_df] = Get_Protein_Correlation_Table(new)
    jo_pro_gene_cor_source.data = jo_Correlation_table_df # update column data source
    jo_pro_cor_data_table.source = jo_pro_gene_cor_source # update table

    kr_pro_gene_cor_source.data = kr_Correlation_table_df  # update column data source
    kr_pro_cor_data_table.source = kr_pro_gene_cor_source  # update table

    me_pro_gene_cor_source.data = me_Correlation_table_df  # update column data source
    me_pro_cor_data_table.source = me_pro_gene_cor_source  # update table

    #protein_table_gene_col_name.title = f'Genes correlated with {new}' # update plot title
    end = time.time()
    print(f"Correlation Table Protein is updated. {end-start}")

def Correlation_Table_mRNA_Update(attrname, old, new):
    start = time.time()
    [jo_m_Correlation_table_df, kr_m_Correlation_table_df, me_m_Correlation_table_df] = Get_mRNA_Correlation_Table(new)
    jo_mrna_gene_cor_source.data = jo_m_Correlation_table_df # update column data source
    jo_mrna_cor_data_table.source = jo_mrna_gene_cor_source # update table

    kr_mrna_gene_cor_source.data = kr_m_Correlation_table_df  # update column data source
    kr_mrna_cor_data_table.source = kr_mrna_gene_cor_source  # update table

    me_mrna_gene_cor_source.data = me_m_Correlation_table_df  # update column data source
    me_mrna_cor_data_table.source = me_mrna_gene_cor_source  # update table

    #mRNA_table_gene_col_name.title = f'Genes correlated with {new}' # update plot title
    end = time.time()
    print(f"Correlation Table mRNA is updated. {end-start}")

################################# Row 1: Protein Complex Subunit Correlation Plot ######################################
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
    #jo_plot_p.circle(subtype_tuple, protein_data, source=johansson_cds[j], color=GENE_COLORS[j], size=4,
    #                 legend_label=GENE_NUMBER[j])
    jo_plot_m.line(subtype_tuple, mRNA_data, source=johansson_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    #jo_plot_m.circle(subtype_tuple, mRNA_data, source=johansson_cds[j], color=GENE_COLORS[j], size=4,
    #                 legend_label=GENE_NUMBER[j])
    TICKER[j] = AutocompleteInput(title=GENE_NUMBER[j], value=INITIAL_GENE[j], width=100, width_policy='fixed',
                                  min_characters=3, completions=Nix(INITIAL_GENE[j], all_unique_genes),
                                  case_sensitive=False)
    TICKER[j].on_change('value', TICKER_FUNCTION_LIST[j])
    TICKER[j].on_change('value', Line_Plot_Update)
    TICKER[j].on_change('value', Scatter_Plot_Update)
    TICKER[j].on_change('value', All_Subtype_Plot_Update)

for j in range(4):
    kr_plot_p.line(subtype_tuple, protein_data, source=krug_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    # kr_plot_p.circle(subtype_tuple, protein_data, source=krug_cds[j], color=GENE_COLORS[j], size=4,
    #                  legend_label=GENE_NUMBER[j])
    kr_plot_m.line(subtype_tuple, mRNA_data, source=krug_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    # kr_plot_m.circle(subtype_tuple, mRNA_data, source=krug_cds[j], color=GENE_COLORS[j], size=4,
    #                  legend_label=GENE_NUMBER[j])
    TICKER[j] = AutocompleteInput(title=GENE_NUMBER[j], value=INITIAL_GENE[j], width=100, width_policy='fixed',
                                  min_characters=3, completions=Nix(INITIAL_GENE[j], all_unique_genes),
                                  case_sensitive=False)
    TICKER[j].on_change('value', TICKER_FUNCTION_LIST[j])
    TICKER[j].on_change('value', Line_Plot_Update)
    TICKER[j].on_change('value', Scatter_Plot_Update)
    TICKER[j].on_change('value', All_Subtype_Plot_Update)

for j in range(4):
    me_plot_p.line(subtype_tuple, protein_data, source=mertins_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    # me_plot_p.circle(subtype_tuple, protein_data, source=mertins_cds[j], color=GENE_COLORS[j], size=4,
    #                  legend_label=GENE_NUMBER[j])
    me_plot_m.line(subtype_tuple, mRNA_data, source=mertins_cds[j], color=GENE_COLORS[j], legend_label=GENE_NUMBER[j])
    # me_plot_m.circle(subtype_tuple, mRNA_data, source=mertins_cds[j], color=GENE_COLORS[j], size=4,
    #                  legend_label=GENE_NUMBER[j])
    TICKER[j] = AutocompleteInput(title=GENE_NUMBER[j], value=INITIAL_GENE[j], width=100, width_policy='fixed',
                                  min_characters=3, completions=Nix(INITIAL_GENE[j], all_unique_genes),
                                  case_sensitive=False)
    TICKER[j].on_change('value', TICKER_FUNCTION_LIST[j])
    TICKER[j].on_change('value', Line_Plot_Update)
    TICKER[j].on_change('value', Scatter_Plot_Update)
    TICKER[j].on_change('value', All_Subtype_Plot_Update)

# styling plots
for i in [jo_plot_p, jo_plot_m, kr_plot_p, kr_plot_m, me_plot_p, me_plot_m]:
    StylePlots(i, PlotID='Correlation')

# Put protein plot and mRNA plot into the same column, with a spacer separating them
line_plot_jo_layout = layout(column(jo_plot_p, Spacer(height=30), jo_plot_m))
line_plot_kr_layout = layout(column(kr_plot_p, Spacer(height=30), kr_plot_m))
line_plot_me_layout = layout(column(me_plot_p, Spacer(height=30), me_plot_m))

# Stack the 3 layouts into tabs
line_plot_tab = Tabs(tabs=[Panel(child=line_plot_jo_layout, title="Johansson"),
                     Panel(child=line_plot_kr_layout, title="Krug"),
                     Panel(child=line_plot_me_layout, title="Mertins")])

#Download data
button1 = Button(label="Download", button_type="success", width=100, width_policy='fixed')  # Gene 1 Data
button1.js_on_event("button_click", CustomJS(args=dict(source=johansson_cds[0]), code=open(join(dirname(__file__),
                                                                                                "Download_Javascript/jo_download.js")).read()))
button1.js_on_event("button_click", CustomJS(args=dict(source=krug_cds[0]), code=open(join(dirname(__file__),
                                                                                           "Download_Javascript/kr_download.js")).read()))
button1.js_on_event("button_click", CustomJS(args=dict(source=mertins_cds[0]), code=open(join(dirname(__file__),
                                                                                              "Download_Javascript/me_download.js")).read()))
####
button2 = Button(label="Download", button_type="success", width=100, width_policy='fixed')  # Gene 2 Data
button2.js_on_event("button_click", CustomJS(args=dict(source=johansson_cds[1]), code=open(join(dirname(__file__),
                                                                                                "Download_Javascript/jo_download.js")).read()))
button2.js_on_event("button_click", CustomJS(args=dict(source=krug_cds[1]), code=open(join(dirname(__file__),
                                                                                           "Download_Javascript/kr_download.js")).read()))
button2.js_on_event("button_click", CustomJS(args=dict(source=mertins_cds[1]), code=open(join(dirname(__file__),
                                                                                              "Download_Javascript/me_download.js")).read()))
####
button3 = Button(label="Download", button_type="success", width=100, width_policy='fixed')  # Gene 3 Data
button3.js_on_event("button_click", CustomJS(args=dict(source=johansson_cds[2]), code=open(join(dirname(__file__),
                                                                                                "Download_Javascript/jo_download.js")).read()))
button3.js_on_event("button_click", CustomJS(args=dict(source=krug_cds[2]), code=open(join(dirname(__file__),
                                                                                           "Download_Javascript/kr_download.js")).read()))
button3.js_on_event("button_click", CustomJS(args=dict(source=mertins_cds[2]), code=open(join(dirname(__file__),
                                                                                              "Download_Javascript/me_download.js")).read()))
####
button4 = Button(label="Download", button_type="success", width=100, width_policy='fixed')  # Gene 4 Data
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

################################# Row 1.5: Pro-Pro & RNA-RNA Correlation Scatter Plot ##################################
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
pprr_scatter_plot_tab = Tabs(tabs=[Panel(child=scatter_plot_jo_layout, title="Johansson"),
                              Panel(child=scatter_plot_kr_layout, title="Krug"),
                              Panel(child=scatter_plot_me_layout, title="Mertins")])

# select widget
x_axis_select = Select(title="X axis", value="Gene 1", options=['Gene 1', 'Gene 2', 'Gene 3', 'Gene 4'], width = 100, width_policy='fixed')
y_axis_select = Select(title="Y axis", value="Gene 2", options=['Gene 1', 'Gene 2', 'Gene 3', 'Gene 4'], width = 100, width_policy='fixed')
selection_update_button = Button(label='Update', button_type='default', width=100, width_policy='fixed')
selection_update_button.on_click(Scatter_Select_Update)

# select widget next to pro-pro/rna-rna scatter plot
select_widget_layout = column(x_axis_select, y_axis_select, selection_update_button)

# static legend under selection dropdown
legend_pic = Div(text="<img src='portal_hdf5/static/Legend pic.png'>")

################################# Row 2: mRNA-Protein Correlation Scatter Plot #########################################
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
scatter_plot_tab = Tabs(tabs=[Panel(child=scatter_plot_jo_layout, title="Johansson"),
                              Panel(child=scatter_plot_kr_layout, title="Krug"),
                              Panel(child=scatter_plot_me_layout, title="Mertins")])

########################################  Row 3: Subtype Mean and SEM plot #############################################
# function call
[jo_subtype_protein_CDS, jo_subtype_mRNA_CDS] = Johansson_Subtype_Avg_Plot_CDS(['ESR1', 'PGR', 'ERBB2', 'MKI67'])
[kr_subtype_protein_CDS, kr_subtype_mRNA_CDS] = Krug_Subtype_Avg_Plot_CDS(['ESR1', 'PGR', 'ERBB2', 'MKI67'])
[me_subtype_protein_CDS, me_subtype_mRNA_CDS] = Mertins_Subtype_Avg_Plot_CDS(['ESR1', 'PGR', 'ERBB2', 'MKI67'])

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
subtype_plot_tab = Tabs(tabs=[Panel(child=subtype_plot_jo_layout, title="Johansson"),
                              Panel(child=subtype_plot_kr_layout, title="Krug"),
                              Panel(child=subtype_plot_me_layout, title="Mertins")])

################################################### Row 4: GxG Correlation #############################################
# function call
[jo_mrna_ERBB2, jo_protein_ERBB2, kr_mrna_ERBB2, kr_protein_ERBB2, me_mrna_ERBB2, me_protein_ERBB2] = Import_Static_Correlation_Table()

# constant
TABLE_WIDTH = 375
TABLE_HEIGHT = 250

# textbox widget and on_change
protein_correlation_textbox = AutocompleteInput(title='Gene:', value=str('ERBB2'), width=180, width_policy='fixed',
                                            min_characters=3, completions = all_unique_genes, case_sensitive=False)
protein_correlation_textbox.on_change('value', Correlation_Table_Protein_Update)
mrna_correlation_textbox = AutocompleteInput(title='Gene:', value=str('ERBB2'), width=180, width_policy='fixed',
                                            min_characters=3, completions = all_unique_genes, case_sensitive=False)
mrna_correlation_textbox.on_change('value', Correlation_Table_mRNA_Update)

# CDS
jo_pro_gene_cor_source = ColumnDataSource(data=jo_protein_ERBB2) # set up columndatasource
kr_pro_gene_cor_source = ColumnDataSource(data=kr_protein_ERBB2)
me_pro_gene_cor_source = ColumnDataSource(data=me_protein_ERBB2)

jo_mrna_gene_cor_source = ColumnDataSource(data=jo_mrna_ERBB2) # set up columndatasource
kr_mrna_gene_cor_source = ColumnDataSource(data=kr_mrna_ERBB2)
me_mrna_gene_cor_source = ColumnDataSource(data=me_mrna_ERBB2)

# table assembly
#protein_table_gene_col_name = TableColumn(field='Gene', title=f'Genes correlated with {protein_correlation_textbox.value}')
#mRNA_table_gene_col_name = TableColumn(field='Gene', title=f'Genes correlated with {mrna_correlation_textbox.value}')

protein_table_columns = [TableColumn(field='Gene', title='Gene'),  # set up mRNA_table_columns
                         TableColumn(field='r', title='Coefficient', formatter=ScientificFormatter(precision=3)),
                         TableColumn(field='p', title='p value', formatter=ScientificFormatter(precision=3))]
mRNA_table_columns = [TableColumn(field='Gene', title='Gene'),  # set up mRNA_table_columns
                      TableColumn(field='r', title='Coefficient', formatter=ScientificFormatter(precision=3)),
                      TableColumn(field='p', title='p value', formatter=ScientificFormatter(precision=3))]

# create table widget by assembling columndatasource and mRNA_table_columns
jo_mrna_cor_data_table = DataTable(source=jo_mrna_gene_cor_source, columns = mRNA_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=False, index_position=None)
kr_mrna_cor_data_table = DataTable(source=kr_mrna_gene_cor_source, columns = mRNA_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=False, index_position=None)
me_mrna_cor_data_table = DataTable(source=me_mrna_gene_cor_source, columns = mRNA_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=False, index_position=None)

jo_pro_cor_data_table = DataTable(source=jo_pro_gene_cor_source, columns = protein_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=False, index_position=None)
kr_pro_cor_data_table = DataTable(source=kr_pro_gene_cor_source, columns = protein_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=False, index_position=None)
me_pro_cor_data_table = DataTable(source=me_pro_gene_cor_source, columns = protein_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=False, index_position=None)

# combine tables into tabs
cor_pro_table_tab = Tabs(tabs=[Panel(child=jo_pro_cor_data_table, title ="Johansson"),
                           Panel(child=kr_pro_cor_data_table, title = "Krug"),
                           Panel(child=me_pro_cor_data_table, title = "Mertins")])

cor_mrna_table_tab = Tabs(tabs=[Panel(child=jo_mrna_cor_data_table, title ="Johansson"),
                           Panel(child=kr_mrna_cor_data_table, title = "Krug"),
                           Panel(child=me_mrna_cor_data_table, title = "Mertins")])

# button and call back
button5 = Button(label="Download", button_type="success", width=180, width_policy='fixed') #mRNA Correlation Table
button5.js_on_event("button_click", CustomJS(args=dict(source=jo_mrna_gene_cor_source),
                                             code=open(join(dirname(__file__),
                                                            "Download_Javascript/jo_cor_table_download.js")).read()))
button5.js_on_event("button_click", CustomJS(args=dict(source=kr_mrna_gene_cor_source),
                                             code=open(join(dirname(__file__),
                                                            "Download_Javascript/kr_cor_table_download.js")).read()))
button5.js_on_event("button_click", CustomJS(args=dict(source=me_mrna_gene_cor_source),
                                             code=open(join(dirname(__file__),
                                                            "Download_Javascript/me_cor_table_download.js")).read()))

button6 = Button(label="Download", button_type="success", width=180, width_policy='fixed') #Protein Correlation Table
button6.js_on_event("button_click", CustomJS(args=dict(source=jo_pro_gene_cor_source),
                                             code=open(join(dirname(__file__),
                                                            "Download_Javascript/jo_cor_table_download.js")).read()))
button6.js_on_event("button_click", CustomJS(args=dict(source=kr_pro_gene_cor_source),
                                             code=open(join(dirname(__file__),
                                                            "Download_Javascript/kr_cor_table_download.js")).read()))
button6.js_on_event("button_click", CustomJS(args=dict(source=me_pro_gene_cor_source),
                                             code=open(join(dirname(__file__),
                                                            "Download_Javascript/me_cor_table_download.js")).read()))

formatted_button5 = row(column(Spacer(height=19),button5))
formatted_button6 = row(column(Spacer(height=19),button6))
########################################################### layout #####################################################
RowSpacer = Spacer(height=30)
PageMargin = Spacer(width=30)

tool_title_text = "Breast Cancer Data Portal"
tool_title_div = Div(text=tool_title_text,
                     style={'font-size': '200%', 'color': 'black', 'font-style': 'normal', 'font-weight': 'bold'},
                     width=1000)

# jo_info_text = '''The Oslo cohort from Johansson et al. has 45 breast cancer tumors including 5 subtypes (basal-like, Her2-enriched, LumA, LumB and normal-like). Each subtype contains 9 samples.
#              The dataset provided identifications of fully quantified 23663 gene transcripts and 9995 proteins.'''
# kr_info_text = '''The BRCA cohort from Krug et al. contains 122 treatment-naive primary breast cancer tumors. Sample subtypes include 29 basal-like, 14 Her2-enriched, 57 LumA, 17 LumB and 5 normal-like samples.
#              The dataset provided identifications of fully quantified 10733 gene transcripts and 7586 proteins.'''
# me_info_text = '''The TCGA-BRCA cohort from Mertins et al. contains 77 breast cancer tumors. Subtypes include 18 basal-like, 12 Her2-enriched, 23 LumA and 24 LumB samples.
#              The dataset provided identifications of fully quantified 14440 gene transcripts, 7016 proteins.'''
# bg_info_div = column(Div(text="Background Information:",
#                          style={'font-size': '120%', 'color': 'black', 'font-style': 'italic', 'font-weight': 'bold'}),
#                      Div(text=jo_info_text,
#                                 style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}),
#                      Div(text=kr_info_text,
#                                 style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}),
#                      Div(text=me_info_text,
#                                 style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}), width=1000)

# line_plot_descriptive_text = "Abundances of proteins that are part of the same complex are tightly correlated across breast tumors. " \
#                              "However, this does not appear to be the case for the corresponding mRNA transcripts. " \
#                              "The plots are initialized showing abundances of structural proteins of complex I of the electron transport chain. " \
#                              "If you wish to hide the curve for a specific gene, click on the corresponding legend entries."
line_plot_div = column(Div(text="Abundance Traces:",
                           style={'font-size': '120%', 'color': 'black', 'font-style': 'italic',
                                       'font-weight': 'bold'}))
                       # Div(text=line_plot_descriptive_text,
                       #     style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}), width=1000)

#scatter_plot_descriptive_text = "Protein and transcript abundances often are not correlated; indicating substantial utilization of post-transcriptional regulatory mechanisms."
scatter_plot_div = column(Div(text="Pairwise Abundances:",
                              style={'font-size': '120%', 'color': 'black', 'font-style': 'italic',
                                      'font-weight': 'bold'}))
                          # Div(text=scatter_plot_descriptive_text,
                          #     style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}), width=1000)

# subtype_plot_text = "Breast cancer subtypes are defined by their gene expression profiles. " \
#                           "The plots are initialized above showing abundances of ER (ESR1), PR (PGR), HER2 (ERBB2), and KI-67 (MKI67); " \
#                           "four immunohistochemical markers commonly used in the clinic. " \
#                           "Values are means +/- standard error of the mean (SEM). " \
#                           "If data are not available for a gene, values will appear as all zeros with no SEMs."
subtype_plot_div = column(Div(text="Abundances by Subtype:",
                              style={'font-size': '120%', 'color': 'black', 'font-style': 'italic', 'font-weight': 'bold'}))
                          # Div(text=subtype_plot_text,
                          #     style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}), width=1000)


protein_cor_table_title_div = column(Div(text="Correlation Table - Protein:",
                                    style={'font-size': '120%', 'color': 'black', 'font-weight': 'bold'}))
mRNA_cor_table_title_div = column(Div(text="Correlation Table - mRNA:",
                                 style={'font-size': '120%', 'color': 'black', 'font-weight': 'bold'}))
# correlation_table_text = "The gene-gene correlation is calculated from protein and mRNA expression data using the pearson correlation. " \
#                "The table shows the top 100 genes that are correlated with the input gene based on the p-value. " \
#                "The computation roughly takes 30 seconds, please be patient."
# correlation_table_div = column(Div(text="Gene-Gene Correlation Table:",
#                                    style={'font-size': '120%', 'color': 'black', 'font-style': 'italic', 'font-weight': 'bold'}),
#                                Div(text=correlation_table_text,
#                                    style={'font-size': '100%', 'color': 'black', 'font-style': 'italic'}), width=1000)

# gene_entry_instruc_text = "To view data for a different gene, type a HGNC gene symbol in the textbox."
# blank_gene_instruc_text = "If you wish to leave certain textboxes blank, please input the corresponding blank entries."
# warning_text = "Please do not enter same gene/blank entry across multiple textboxes."
# download_text = 'Click "download" to save the protein and mRNA data for corresponding genes from the three studies. ' \
#               'Please use Firefox or Google Chrome for better download experience.'
# red_text_div = column(Div(text=gene_entry_instruc_text, style={'font-size': '100%', 'color': 'red', 'font-style': 'italic'}, width=300),
#                     Div(text=blank_gene_instruc_text, style={'font-size': '100%', 'color': 'red', 'font-style': 'italic'}, width=300),
#                     Div(text=warning_text, style={'font-size': '100%', 'color': 'red', 'font-style': 'italic'}, width=300),
#                     Div(text=download_text, style={'font-size': '100%', 'color': 'red', 'font-style': 'italic'}, width=300))

gene_hide_text_div = Div(text="If you wish to hide certain genes, please click on the interactive legend.",
                         style={'font-size': '100%', 'color': 'red'}, width=200, width_policy='fixed')

gene_entry_section = column(tickers_buttons_layout)
line_plot_section = column(line_plot_div, row(line_plot_tab, column(gene_entry_section, gene_hide_text_div)))
select_widget_section = column(select_widget_layout)
scatter_plot_section = column(scatter_plot_div, row(scatter_plot_tab, Spacer(width=60), pprr_scatter_plot_tab, column(select_widget_section, legend_pic)))
subtype_plot_section = row(column(subtype_plot_div, subtype_plot_tab),column(Spacer(height=60),legend_pic))
correlation_table_section = column(row(column(protein_cor_table_title_div, row(protein_correlation_textbox, formatted_button6), cor_pro_table_tab),
                          column(mRNA_cor_table_title_div, row(mrna_correlation_textbox,formatted_button5), cor_mrna_table_tab)))

l = layout([
    [PageMargin, tool_title_div],
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
print(f'total response time: {elapsed}')

