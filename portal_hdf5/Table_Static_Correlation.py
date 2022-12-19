# initiate correlation table (static)

from bokeh.models import ColumnDataSource, TableColumn, DataTable, ScientificFormatter
from Import_Files import Import_Static_Correlation_Table

# function call, import the static correlation tables for ERBB2
[jo_mrna_ERBB2, jo_protein_ERBB2, kr_mrna_ERBB2, kr_protein_ERBB2, me_mrna_ERBB2, me_protein_ERBB2, summary_protein_ERBB2, summary_mrna_ERBB2] = Import_Static_Correlation_Table()

# constants for datatable widget width and height
TABLE_WIDTH = 900
TABLE_HEIGHT = 250

# static tables ColumnDataSource set up
jo_pro_gene_cor_source = ColumnDataSource(data=jo_protein_ERBB2) # set up columndatasource
kr_pro_gene_cor_source = ColumnDataSource(data=kr_protein_ERBB2)
me_pro_gene_cor_source = ColumnDataSource(data=me_protein_ERBB2)
summary_pro_source = ColumnDataSource(data=summary_protein_ERBB2)

jo_mrna_gene_cor_source = ColumnDataSource(data=jo_mrna_ERBB2) # set up columndatasource
kr_mrna_gene_cor_source = ColumnDataSource(data=kr_mrna_ERBB2)
me_mrna_gene_cor_source = ColumnDataSource(data=me_mrna_ERBB2)
summary_mrna_source = ColumnDataSource(data=summary_mrna_ERBB2)

# setting up column names, this is the version of not updating the gene column name
protein_table_gene_col_name = TableColumn(field='Gene', title='Gene')
mRNA_table_gene_col_name = TableColumn(field='Gene', title='Gene')

protein_table_columns = [protein_table_gene_col_name,  # set up protein_table_columns
                         TableColumn(field='r', title='Coefficient', formatter=ScientificFormatter(precision=3)),
                         TableColumn(field='p', title='p value', formatter=ScientificFormatter(precision=3)),
                         TableColumn(field='Name', title='Name')]
mRNA_table_columns = [mRNA_table_gene_col_name,  # set up mRNA_table_columns
                      TableColumn(field='r', title='Coefficient', formatter=ScientificFormatter(precision=3)),
                      TableColumn(field='p', title='p value', formatter=ScientificFormatter(precision=3)),
                      TableColumn(field='Name', title='Name')]

protein_summary_table_columns = [protein_table_gene_col_name, # set up protein_summary_table_columns
                         TableColumn(field='p avg', title='Averaged p values'),
                         TableColumn(field='Name', title='Name')]

mRNA_summary_table_columns = [mRNA_table_gene_col_name, # set up mRNA_summary_table_columns
                         TableColumn(field='p avg', title='Averaged p values'),
                         TableColumn(field='Name', title='Name')]

# this function is unknown, probably not useful
def Correlation_DataTable():
    jo_pro_cor_data_table = DataTable(source=jo_pro_gene_cor_source, columns=protein_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, editable=True, index_position=None)
    kr_pro_cor_data_table = DataTable(source=kr_pro_gene_cor_source, columns=protein_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, editable=True, index_position=None)
    me_pro_cor_data_table = DataTable(source=me_pro_gene_cor_source, columns=protein_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, editable=True, index_position=None)
    summary_pro_data_table = DataTable(source=summary_pro_source, columns=protein_summary_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, editable=True, index_position=None)

    jo_mrna_cor_data_table = DataTable(source=jo_mrna_gene_cor_source, columns=mRNA_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, editable=True, index_position=None)
    kr_mrna_cor_data_table = DataTable(source=kr_mrna_gene_cor_source, columns=mRNA_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, editable=True, index_position=None)
    me_mrna_cor_data_table = DataTable(source=me_mrna_gene_cor_source, columns=mRNA_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, editable=True, index_position=None)
    summary_mrna_data_table = DataTable(source=summary_mrna_source, columns=mRNA_summary_table_columns, width=TABLE_WIDTH, height=TABLE_HEIGHT, editable=True, index_position=None)

    return jo_pro_cor_data_table, kr_pro_cor_data_table, me_pro_cor_data_table, summary_pro_data_table,\
           jo_mrna_cor_data_table, kr_mrna_cor_data_table, me_mrna_cor_data_table, summary_mrna_data_table