############################################ Import ####################################################################
from bokeh.io import show
from bokeh.layouts import layout, column, row
from bokeh.models import FactorRange, Spacer, Tabs, Panel, AutocompleteInput, ColumnDataSource, TableColumn, DataTable
from bokeh.plotting import figure
from bokeh.transform import factor_cmap

from Plot_All_Styling import StylePlots
from Line_Plot import Johansson_Line_Plot, Krug_Line_Plot, Mertins_Line_Plot
from Scatter_Plot import Pro_Pro_Scatter_Plot, RNA_RNA_Scatter_Plot, Johansson_RNA_Pro_Scatter_Plot, Krug_RNA_Pro_Scatter_Plot, Mertins_RNA_Pro_Scatter_Plot
from Subtype_Plot import Johansson_Subtype_Plot, Krug_Subtype_Plot, Mertins_Subtype_Plot
from Table_Correlation import Correlation_DataTable

import time

################################# Row 1: Protein Complex Subunit Correlation Plot ######################################
start_time = time.time()

# function call
[jo_plot_p, jo_plot_m] = Johansson_Line_Plot()
[kr_plot_p, kr_plot_m] = Krug_Line_Plot()
[me_plot_p, me_plot_m] = Mertins_Line_Plot()

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

################################# Row 1.5: Pro-Pro & RNA-RNA Correlation Scatter Plot ##################################
# function call to get the plots
[jo_pro_pro_scatter_plot, kr_pro_pro_scatter_plot, me_pro_pro_scatter_plot] = Pro_Pro_Scatter_Plot()
[jo_rna_rna_scatter_plot, kr_rna_rna_scatter_plot, me_rna_rna_scatter_plot] = RNA_RNA_Scatter_Plot()

# put plots into a list and run through styling function
p_p_plots = [jo_pro_pro_scatter_plot, kr_pro_pro_scatter_plot, me_pro_pro_scatter_plot]
r_r_plots = [jo_rna_rna_scatter_plot, kr_rna_rna_scatter_plot, me_rna_rna_scatter_plot]
for i,j in zip(p_p_plots,r_r_plots):
    StylePlots(i, PlotID='mRNA-Prot')
    StylePlots(j, PlotID='mRNA-Prot')

# layout
scatter_plot_jo_layout = layout(column(row(p_p_plots[0],r_r_plots[0])))

scatter_plot_kr_layout = layout(column(row(p_p_plots[1],r_r_plots[1])))

scatter_plot_me_layout = layout(column(row(p_p_plots[2],r_r_plots[2])))

# stack 3 layouts into tab
pprr_scatter_plot_tab = Tabs(tabs=[Panel(child=scatter_plot_jo_layout, title="Johansson"),
                              Panel(child=scatter_plot_kr_layout, title="Krug"),
                              Panel(child=scatter_plot_me_layout, title="Mertins")])

################################# Row 2: mRNA-Protein Correlation Scatter Plot #########################################
# function calls
jo_plot_mRNA_prot = Johansson_RNA_Pro_Scatter_Plot()
kr_plot_mRNA_prot = Krug_RNA_Pro_Scatter_Plot()
me_plot_mRNA_prot = Mertins_RNA_Pro_Scatter_Plot()

# stylish plots
for i in jo_plot_mRNA_prot:
    StylePlots(i, PlotID='mRNA-Prot')
for i in kr_plot_mRNA_prot:
    StylePlots(i, PlotID='mRNA-Prot')
for i in me_plot_mRNA_prot:
    StylePlots(i, PlotID='mRNA-Prot')

# put gene1 and gene2 mRNA-protein correlation plots in the same row, with a spacer between
# put gene3 and gene4 mRNA-protein correlation plots in the row below, with a spacer between.
# spacer to divide the rows
scatter_plot_jo_layout = layout(column(row(jo_plot_mRNA_prot[0],jo_plot_mRNA_prot[1]),
                                       Spacer(height=30),
                                       row(jo_plot_mRNA_prot[2],jo_plot_mRNA_prot[3])))
scatter_plot_kr_layout = layout(column(row(kr_plot_mRNA_prot[0],kr_plot_mRNA_prot[1]),
                                       Spacer(height=30),
                                       row(kr_plot_mRNA_prot[2],kr_plot_mRNA_prot[3])))
scatter_plot_me_layout = layout(column(row(me_plot_mRNA_prot[0],me_plot_mRNA_prot[1]),
                                       Spacer(height=30),
                                       row(me_plot_mRNA_prot[2],me_plot_mRNA_prot[3])))
# stack 3 layouts into tab
scatter_plot_tab = Tabs(tabs=[Panel(child=scatter_plot_jo_layout, title="Johansson"),
                              Panel(child=scatter_plot_kr_layout, title="Krug"),
                              Panel(child=scatter_plot_me_layout, title="Mertins")])

########################################  Row 3: Subtype Mean and SEM plot #############################################
# function call
[jo_protein_subtype_plot, jo_mRNA_subtype_plot] = Johansson_Subtype_Plot()
[kr_protein_subtype_plot, kr_mRNA_subtype_plot] = Krug_Subtype_Plot()
[me_protein_subtype_plot, me_mRNA_subtype_plot] = Mertins_Subtype_Plot()

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
[jo_pro_cor_data_table, kr_pro_cor_data_table, me_pro_cor_data_table, jo_mrna_cor_data_table, kr_mrna_cor_data_table, me_mrna_cor_data_table] = Correlation_DataTable()

# textbox widget and on_change
#mrna_correlation_textbox = AutocompleteInput(title='Gene Name for mRNA correlation table:', value=str('ERBB2'), width=150, width_policy ='auto',
#                                             min_characters=3, completions = mrna_cor_gene_list, case_sensitive=False)
#mrna_correlation_textbox.on_change('value', update_mrna_correlation)
#
#protein_correlation_textbox = AutocompleteInput(title='Gene Name for protein correlation table:', value=str('ERBB2'), width=150, width_policy ='auto',
#                                             min_characters=3, completions = protein_cor_gene_list, case_sensitive=False)
#protein_correlation_textbox.on_change('value', update_protein_correlation)


# combine tables into tabs
cor_pro_table_tab = Tabs(tabs=[Panel(child=jo_pro_cor_data_table, title ="Johansson"),
                           Panel(child=kr_pro_cor_data_table, title = "Krug"),
                           Panel(child=me_pro_cor_data_table, title = "Mertins")])

cor_mrna_table_tab = Tabs(tabs=[Panel(child=jo_mrna_cor_data_table, title ="Johansson"),
                           Panel(child=kr_mrna_cor_data_table, title = "Krug"),
                           Panel(child=me_mrna_cor_data_table, title = "Mertins")])

# show
l=layout(line_plot_tab, pprr_scatter_plot_tab, scatter_plot_tab, subtype_plot_tab, row(cor_pro_table_tab, cor_mrna_table_tab))
show(l)
end_time=time.time()
elapsed=end_time-start_time
print(f'total time: {elapsed}')

