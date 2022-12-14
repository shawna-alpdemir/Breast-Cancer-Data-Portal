# function that creates subtype plot CDS

#from pdb import set_trace
from bokeh.models import ColumnDataSource
from portal_hdf5.Archive.Subtype_Average_DF import Johansson_Subtype_Avg_SEM_DFs, Krug_Subtype_Avg_SEM_DFs, Mertins_Subtype_Avg_SEM_DFs


########################################################################################################################
##################################################### Johansson ########################################################

def Johansson_Subtype_Avg_Plot_CDS(genelist):
    """ Create ColumnDataSource dictionary for plotting subtype mean and SEM purposes """
    # Create empty list to store data iterated from the for loop below
    x_axis_subtype = []

    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    subtype = []

    for gene in genelist:
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

        subtype.extend(['Basal', 'Her2', 'LumA', 'LumB', 'Norm'])
    # create the plot data sources from the dictionaries for each gene.
    jo_source_dict_subtype_protein = {'x': x_axis_subtype, 'y': y_axis_subtype_protein,
              'upper': upper_bar_protein, 'lower': lower_bar_protein, 'legend_group': subtype}

    jo_subtype_protein_CDS = ColumnDataSource(data=jo_source_dict_subtype_protein)

    jo_source_dict_subtype_mrna = {'x': x_axis_subtype, 'y': y_axis_subtype_mrna,
              'upper': upper_bar_mrna, 'lower': lower_bar_mrna, 'legend_group': subtype}

    jo_subtype_mRNA_CDS = ColumnDataSource(data=jo_source_dict_subtype_mrna)

    # The data dictionary will look like this:
    # 'x': [(geneA, basal), (geneA, her2)....]
    # 'y': [mean for geneA x basal, mean for geneA x her2...]
    # 'upper' [mean+SEM for geneA x basal, mean+SEM for geneA x her2...]
    # 'lower' [mean-SEM for geneA x basal, mean-SEM for geneA x her2...]
    # 'legend_group = [5 subtypes * 4]

    return jo_subtype_protein_CDS, jo_subtype_mRNA_CDS


########################################################################################################################
##################################################### Krug #############################################################

def Krug_Subtype_Avg_Plot_CDS(genelist):
    """ Create ColumnDataSource dictionary for plotting subtype mean and SEM purposes """
    # Create empty list to store data iterated from the for loop below
    x_axis_subtype = []

    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    subtype = []

    for gene in genelist:
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

        subtype.extend(['Basal', 'Her2', 'LumA', 'LumB', 'Norm'])
    # create the plot data sources from the dictionaries for each gene.
    kr_source_dict_subtype_protein = {'x': x_axis_subtype, 'y': y_axis_subtype_protein,
              'upper': upper_bar_protein, 'lower': lower_bar_protein, 'legend_group': subtype}

    kr_subtype_protein_CDS = ColumnDataSource(data=kr_source_dict_subtype_protein)

    kr_source_dict_subtype_mrna = {'x': x_axis_subtype, 'y': y_axis_subtype_mrna,
              'upper': upper_bar_mrna, 'lower': lower_bar_mrna, 'legend_group': subtype}

    kr_subtype_mRNA_CDS = ColumnDataSource(data=kr_source_dict_subtype_mrna)

    # The data dictionary will look like this:
    # 'x': [(geneA, basal), (geneA, her2)....]
    # 'y': [mean for geneA x basal, mean for geneA x her2...]
    # 'upper' [mean+SEM for geneA x basal, mean+SEM for geneA x her2...]
    # 'lower' [mean-SEM for geneA x basal, mean-SEM for geneA x her2...]
    # 'legend_group = [5 subtypes * 4]

    return kr_subtype_protein_CDS, kr_subtype_mRNA_CDS

########################################################################################################################
##################################################### Mertins ##########################################################

def Mertins_Subtype_Avg_Plot_CDS(genelist):
    """ Create ColumnDataSource dictionary for plotting subtype mean and SEM purposes """
    # Create empty list to store data iterated from the for loop below
    x_axis_subtype = []

    y_axis_subtype_protein = []
    upper_bar_protein = []
    lower_bar_protein = []

    y_axis_subtype_mrna = []
    upper_bar_mrna = []
    lower_bar_mrna = []

    subtype = []

    for gene in genelist:
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

        subtype.extend(['Basal', 'Her2', 'LumA', 'LumB'])
    # create the plot data sources from the dictionaries for each gene.
    me_source_dict_subtype_protein = {'x': x_axis_subtype, 'y': y_axis_subtype_protein,
              'upper': upper_bar_protein, 'lower': lower_bar_protein, 'legend_group': subtype}

    me_subtype_protein_CDS = ColumnDataSource(data=me_source_dict_subtype_protein)

    me_source_dict_subtype_mrna = {'x': x_axis_subtype, 'y': y_axis_subtype_mrna,
              'upper': upper_bar_mrna, 'lower': lower_bar_mrna, 'legend_group': subtype}

    me_subtype_mRNA_CDS = ColumnDataSource(data=me_source_dict_subtype_mrna)

    # The data dictionary will look like this:
    # 'x': [(geneA, basal), (geneA, her2)....]
    # 'y': [mean for geneA x basal, mean for geneA x her2...]
    # 'upper' [mean+SEM for geneA x basal, mean+SEM for geneA x her2...]
    # 'lower' [mean-SEM for geneA x basal, mean-SEM for geneA x her2...]
    # 'legend_group = [5 subtypes * 4]

    return me_subtype_protein_CDS, me_subtype_mRNA_CDS

