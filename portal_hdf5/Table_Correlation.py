from bokeh.models import ColumnDataSource, TableColumn, DataTable
from Import_Files import Import_Static_Correlation_Table

# function call
[jo_mrna_ERBB2, jo_protein_ERBB2, kr_mrna_ERBB2, kr_protein_ERBB2, me_mrna_ERBB2, me_protein_ERBB2] = Import_Static_Correlation_Table()

# constant
TABLE_WIDTH = 500
TABLE_HEIGHT = 300

# aesthetic --  static tables
jo_pro_gene_cor_source = ColumnDataSource(data=jo_protein_ERBB2) # set up columndatasource
kr_pro_gene_cor_source = ColumnDataSource(data=kr_protein_ERBB2)
me_pro_gene_cor_source = ColumnDataSource(data=me_protein_ERBB2)

jo_mrna_gene_cor_source = ColumnDataSource(data=jo_mrna_ERBB2) # set up columndatasource
kr_mrna_gene_cor_source = ColumnDataSource(data=kr_mrna_ERBB2)
me_mrna_gene_cor_source = ColumnDataSource(data=me_mrna_ERBB2)

protein_table_gene_col_name = TableColumn(field='Gene', title=f'Genes correlated with ERBB2')
mRNA_table_gene_col_name = TableColumn(field='Gene', title=f'Genes correlated with ERBB2')

protein_table_columns = [protein_table_gene_col_name,  # set up mRNA_table_columns
                         TableColumn(field='r', title='Coefficient'),
                         TableColumn(field='p', title='p value')]
mRNA_table_columns = [mRNA_table_gene_col_name,  # set up mRNA_table_columns
                      TableColumn(field='r', title='Coefficient'),
                      TableColumn(field='p', title='p value')]

# create table widget by assembling columndatasource and mRNA_table_columns
def Correlation_DataTable():
    jo_mrna_cor_data_table = DataTable(source=jo_mrna_gene_cor_source, columns = mRNA_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)
    kr_mrna_cor_data_table = DataTable(source=kr_mrna_gene_cor_source, columns = mRNA_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)
    me_mrna_cor_data_table = DataTable(source=me_mrna_gene_cor_source, columns = mRNA_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)

    jo_pro_cor_data_table = DataTable(source=jo_pro_gene_cor_source, columns = protein_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)
    kr_pro_cor_data_table = DataTable(source=kr_pro_gene_cor_source, columns = protein_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)
    me_pro_cor_data_table = DataTable(source=me_pro_gene_cor_source, columns = protein_table_columns, width=TABLE_WIDTH, height = TABLE_HEIGHT, editable=True, index_position=None)

    return jo_pro_cor_data_table, kr_pro_cor_data_table, me_pro_cor_data_table, jo_mrna_cor_data_table, kr_mrna_cor_data_table, me_mrna_cor_data_table