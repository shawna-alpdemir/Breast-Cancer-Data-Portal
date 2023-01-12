# function that creates subtype plot CDS
from pdb import set_trace
from bokeh.models import ColumnDataSource

def Johansson_Subtype_Avg_Plot_CDS():
    jo_source_dict_subtype_protein = {'x': [('ESR1', 'Basal'), ('ESR1', 'Her2'), ('ESR1', 'LumA'), ('ESR1', 'LumB'), ('ESR1', 'Norm'), ('PGR', 'Basal'), ('PGR', 'Her2'), ('PGR', 'LumA'), ('PGR', 'LumB'), ('PGR', 'Norm'), ('ERBB2', 'Basal'), ('ERBB2', 'Her2'), ('ERBB2', 'LumA'), ('ERBB2', 'LumB'), ('ERBB2', 'Norm'), ('MKI67', 'Basal'), ('MKI67', 'Her2'), ('MKI67', 'LumA'), ('MKI67', 'LumB'), ('MKI67', 'Norm')],
                                    'y': [-0.973679125, -0.545407126, 0.857246294, 0.829601229, -0.167761273, -0.581799977, -0.287072222, 0.52324874, 0.461697734, -0.116074275, -0.57581931, 1.521765014, -0.308505996, -0.319949889, -0.31748982, 0.93238547, 0.020248225, -0.636275203, 0.111467467, -0.427825959],
                                    'upper': [-0.929665278, -0.471794554, 1.110458656, 1.175183177, 0.030791655, -0.462278129, -0.142204265, 0.886350399, 0.905768487, 0.092619902, -0.528359393, 1.964042697, -0.272766583, -0.290251543, -0.243932496, 1.4781631800000001, 0.17227413, -0.5695743839999999, 0.270682026, -0.380808343],
                                    'lower': [-1.017692972, -0.619019698, 0.6040339319999999, 0.48401928099999997, -0.366314201, -0.701321825, -0.431940179, 0.16014708099999997, 0.01762698100000004, -0.32476845200000004, -0.6232792269999999, 1.0794873310000002, -0.344245409, -0.349648235, -0.39104714399999996, 0.3866077600000001, -0.13177767999999998, -0.702976022, -0.047747092000000005, -0.47484357499999996],
                                    'legend_group': ['Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm']}
    jo_subtype_protein_CDS = ColumnDataSource(data=jo_source_dict_subtype_protein)

    jo_source_dict_subtype_mrna = {'x': [('ESR1', 'Basal'), ('ESR1', 'Her2'), ('ESR1', 'LumA'), ('ESR1', 'LumB'), ('ESR1', 'Norm'), ('PGR', 'Basal'), ('PGR', 'Her2'), ('PGR', 'LumA'), ('PGR', 'LumB'), ('PGR', 'Norm'), ('ERBB2', 'Basal'), ('ERBB2', 'Her2'), ('ERBB2', 'LumA'), ('ERBB2', 'LumB'), ('ERBB2', 'Norm'), ('MKI67', 'Basal'), ('MKI67', 'Her2'), ('MKI67', 'LumA'), ('MKI67', 'LumB'), ('MKI67', 'Norm')],
                               'y': [-1.202715768, -0.434802393, 0.899017477, 0.8198445, -0.081343817, -1.110246481, -0.363265768, 0.671865036, 0.561091651, 0.240555563, -0.229170399, 1.363879049, -0.532954008, -0.605044747, 0.003290105, 1.056041101, 0.492141645, -1.029462843, 0.268199034, -0.786918937],
                               'upper': [-1.08250957, -0.19674832, 1.025985096, 0.906391559, 0.20484002, -1.057584071, -0.086750206, 0.882136493, 0.848133953, 0.49388748699999996, -0.07838340599999999, 1.748470206, -0.45772243500000004, -0.547470634, 0.233476416, 1.259419752, 0.690937985, -0.8308254179999999, 0.393658479, -0.5820533210000001],
                               'lower': [-1.322921966, -0.672856466, 0.772049858, 0.7332974409999999, -0.367527654, -1.1629088909999998, -0.6397813299999999, 0.46159357900000003, 0.274049349, -0.012776360999999986, -0.379957392, 0.9792878919999999, -0.6081855810000001, -0.66261886, -0.22689620600000002, 0.85266245, 0.293345305, -1.228100268, 0.142739589, -0.991784553],
                               'legend_group': ['Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm']}
    jo_subtype_mRNA_CDS = ColumnDataSource(data=jo_source_dict_subtype_mrna)

    return jo_subtype_protein_CDS, jo_subtype_mRNA_CDS

def Krug_Subtype_Avg_Plot_CDS():
    kr_source_dict_subtype_protein = {'x': [('ESR1', 'Basal'), ('ESR1', 'Her2'), ('ESR1', 'LumA'), ('ESR1', 'LumB'), ('ESR1', 'Norm'), ('PGR', 'Basal'), ('PGR', 'Her2'), ('PGR', 'LumA'), ('PGR', 'LumB'), ('PGR', 'Norm'), ('ERBB2', 'Basal'), ('ERBB2', 'Her2'), ('ERBB2', 'LumA'), ('ERBB2', 'LumB'), ('ERBB2', 'Norm'), ('MKI67', 'Basal'), ('MKI67', 'Her2'), ('MKI67', 'LumA'), ('MKI67', 'LumB'), ('MKI67', 'Norm')],
                                  'y': [-0.881337298, -0.939449733, 0.51435491, 0.819173402, -0.906619961, -0.661662884, -0.652878504, 0.503549828, 0.223666976, -0.835231222, -0.259523372, 1.254511038, -0.186337085, -0.178408828, 0.723437435, 1.234126928, 0.035030251, -0.649001676, 0.221208868, -0.609511931],
                                  'upper': [-0.796351692, -0.744872287, 0.620882665, 0.982918007, -0.795029452, -0.570489781, -0.438770903, 0.623791667, 0.514774103, -0.723136376, -0.159547691, 1.691158415, -0.096592627, 0.102033131, 1.108110868, 1.393035503, 0.207591995, -0.576455772, 0.385985119, -0.385077082],
                                  'lower': [-0.966322904, -1.134027179, 0.407827155, 0.655428797, -1.01821047, -0.752835987, -0.866986105, 0.383307989, -0.067440151, -0.947326068, -0.359499053, 0.817863661, -0.276081543, -0.458850787, 0.338764002, 1.075218353, -0.137531493, -0.72154758, 0.056432617, -0.83394678],
                                  'legend_group': ['Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm']}
    kr_subtype_protein_CDS = ColumnDataSource(data=kr_source_dict_subtype_protein)

    kr_source_dict_subtype_mrna = {'x': [('ESR1', 'Basal'), ('ESR1', 'Her2'), ('ESR1', 'LumA'), ('ESR1', 'LumB'), ('ESR1', 'Norm'), ('PGR', 'Basal'), ('PGR', 'Her2'), ('PGR', 'LumA'), ('PGR', 'LumB'), ('PGR', 'Norm'), ('ERBB2', 'Basal'), ('ERBB2', 'Her2'), ('ERBB2', 'LumA'), ('ERBB2', 'LumB'), ('ERBB2', 'Norm'), ('MKI67', 'Basal'), ('MKI67', 'Her2'), ('MKI67', 'LumA'), ('MKI67', 'LumB'), ('MKI67', 'Norm')],
                               'y': [-1.241516223, -1.043948988, 0.687274447, 0.800030553, -0.431181309, -1.171341538, -0.737404845, 0.736663809, 0.220556472, -0.289344944, -0.688779177, 1.182127128, 0.023602744, -0.035783011, 0.537554223, 0.902412891, 0.209631299, -0.64044023, 0.503312448, -0.231206108],
                               'upper': [-1.147046663, -0.878074833, 0.73277988, 0.8477071, -0.055258773, -1.081769362, -0.601894303, 0.812692502, 0.432752517, -0.081388514, -0.578924678, 1.604106075, 0.108180444, 0.204858358, 1.056576661, 1.034871707, 0.357468474, -0.522738019, 0.609854484, 0.173890693],
                               'lower': [-1.335985783, -1.209823143, 0.641769014, 0.752354006, -0.807103845, -1.260913714, -0.872915387, 0.660635116, 0.008360427, -0.497301374, -0.798633676, 0.760148181, -0.060974956, -0.27642438, 0.018531785, 0.769954075, 0.061794124, -0.758142441, 0.396770412, -0.636302909],
                               'legend_group': ['Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm', 'Basal', 'Her2', 'LumA', 'LumB', 'Norm']}
    kr_subtype_mRNA_CDS = ColumnDataSource(data=kr_source_dict_subtype_mrna)

    return kr_subtype_protein_CDS, kr_subtype_mRNA_CDS

def Mertins_Subtype_Avg_Plot_CDS():
    me_source_dict_subtype_protein = {'x': [('ESR1', 'Basal'), ('ESR1', 'Her2'), ('ESR1', 'LumA'), ('ESR1', 'LumB'), ('PGR', 'Basal'), ('PGR', 'Her2'), ('PGR', 'LumA'), ('PGR', 'LumB'), ('ERBB2', 'Basal'), ('ERBB2', 'Her2'), ('ERBB2', 'LumA'), ('ERBB2', 'LumB'), ('MKI67', 'Basal'), ('MKI67', 'Her2'), ('MKI67', 'LumA'), ('MKI67', 'LumB')],
                                   'y': [-0.891419303, -0.752694756, 0.403585818, 0.658142113, -0.646942425, -0.619333737, 0.5581647, 0.25996585, -0.416929282, 1.208160657, -0.359746797, 0.05337398, 0.74359176, 0.033204967, -0.693254032, 0.090072143],
                                   'upper': [-0.747880267, -0.603403796, 0.566650761, 0.843937126, -0.509468596, -0.425742751, 0.773510396, 0.459060951, -0.295528327, 1.559479085, -0.263533628, 0.27788749, 0.991365353, 0.272670476, -0.524425662, 0.254959774],
                                   'lower': [-1.034958339, -0.901985716, 0.240520875, 0.4723471, -0.784416254, -0.812924723, 0.342819004, 0.060870749, -0.538330237, 0.856842229, -0.455959966, -0.17113953, 0.495818167, -0.206260542, -0.862082402, -0.074815488],
                                   'legend_group': ['Basal', 'Her2', 'LumA', 'LumB', 'Basal', 'Her2', 'LumA', 'LumB', 'Basal', 'Her2', 'LumA', 'LumB', 'Basal', 'Her2', 'LumA', 'LumB']}
    me_subtype_protein_CDS = ColumnDataSource(data=me_source_dict_subtype_protein)

    me_source_dict_subtype_mrna = {'x': [('ESR1', 'Basal'), ('ESR1', 'Her2'), ('ESR1', 'LumA'), ('ESR1', 'LumB'), ('PGR', 'Basal'), ('PGR', 'Her2'), ('PGR', 'LumA'), ('PGR', 'LumB'), ('ERBB2', 'Basal'), ('ERBB2', 'Her2'), ('ERBB2', 'LumA'), ('ERBB2', 'LumB'), ('MKI67', 'Basal'), ('MKI67', 'Her2'), ('MKI67', 'LumA'), ('MKI67', 'LumB')],
                               'y': [-1.315851054, -0.702931457, 0.631147787, 0.733504056, -0.93308612, -0.750719542, 0.852655218, 0.258046443, -0.665570021, 1.448408904, -0.141131598, -0.089775821, 0.402392073, 0.521426459, -0.655263902, 0.065453956],
                               'upper': [-1.182581908, -0.501495987, 0.69924904, 0.815287628, -0.803402169, -0.561682036, 0.960470081, 0.446001972, -0.553339946, 1.758251228, -0.068049279, 0.114219506, 0.670109184, 0.709700436, -0.473532407, 0.238462098],
                               'lower': [-1.4491202, -0.904366927, 0.563046534, 0.651720484, -1.062770071, -0.939757048, 0.744840355, 0.070090914, -0.777800096, 1.13856658, -0.214213917, -0.293771148, 0.134674962, 0.333152482, -0.836995397, -0.107554186],
                               'legend_group': ['Basal', 'Her2', 'LumA', 'LumB', 'Basal', 'Her2', 'LumA', 'LumB', 'Basal', 'Her2', 'LumA', 'LumB', 'Basal', 'Her2', 'LumA', 'LumB']}
    me_subtype_mRNA_CDS = ColumnDataSource(data=me_source_dict_subtype_mrna)

    return me_subtype_protein_CDS, me_subtype_mRNA_CDS

# methods to generate CDS

# def Subtype_Avg_Plot_CDS(genelist):
#     """ Create ColumnDataSource dictionary for plotting subtype mean and SEM purposes """
#     # Create empty list to store data iterated from the for loop below
#     x_axis_subtype = []
#
#     y_axis_subtype_protein = []
#     upper_bar_protein = []
#     lower_bar_protein = []
#
#     y_axis_subtype_mrna = []
#     upper_bar_mrna = []
#     lower_bar_mrna = []
#
#     subtype = []
#
#     for gene in genelist:
#         try:
#             gene_sub_avg_pro = me_pro.loc[gene,['basal_avg','her2_avg','lumA_avg','lumB_avg']].tolist()
#             gene_sub_upper_pro = me_pro.loc[gene,['basal_up','her2_up','lumA_up','lumB_up']].tolist()
#             gene_sub_lower_pro = me_pro.loc[gene,['basal_down','her2_down','lumA_down','lumB_down']].tolist()
#
#             gene_sub_avg_rna = me_rna.loc[gene,['basal_avg','her2_avg','lumA_avg','lumB_avg']].tolist()
#             gene_sub_upper_rna = me_rna.loc[gene,['basal_up','her2_up','lumA_up','lumB_up']].tolist()
#             gene_sub_lower_rna = me_rna.loc[gene,['basal_down','her2_down','lumA_down','lumB_down']].tolist()
#
#         except ValueError:
#             gene_sub_avg_pro = [0, 0, 0, 0, 0]
#             gene_sub_upper_pro = [0, 0, 0, 0, 0]
#             gene_sub_lower_pro = [0, 0, 0, 0, 0]
#             gene_sub_avg_rna = [0, 0, 0, 0, 0]
#             gene_sub_upper_rna = [0, 0, 0, 0, 0]
#             gene_sub_lower_rna = [0, 0, 0, 0, 0]
#
#         x_axis_subtype.extend([(gene,x) for x in ['Basal', 'Her2', 'LumA', 'LumB']])
#
#         y_axis_subtype_protein.extend(gene_sub_avg_pro)
#         y_axis_subtype_mrna.extend(gene_sub_avg_rna)
#
#         upper_bar_protein.extend(float(x) for x in gene_sub_upper_pro)
#         lower_bar_protein.extend(float(x) for x in gene_sub_lower_pro)
#
#         upper_bar_mrna.extend(float(x) for x in gene_sub_upper_rna)
#         lower_bar_mrna.extend(float(x) for x in gene_sub_lower_rna)
#
#         subtype.extend(['Basal', 'Her2', 'LumA', 'LumB'])
#
#     # create the plot data sources from the dictionaries for each gene.
#     source_dict_subtype_protein = {'x': x_axis_subtype, 'y': y_axis_subtype_protein,
#               'upper': upper_bar_protein, 'lower': lower_bar_protein, 'legend_group': subtype}
#
#     subtype_protein_CDS = ColumnDataSource(data=source_dict_subtype_protein)
#
#     source_dict_subtype_mrna = {'x': x_axis_subtype, 'y': y_axis_subtype_mrna,
#               'upper': upper_bar_mrna, 'lower': lower_bar_mrna, 'legend_group': subtype}
#
#     subtype_mRNA_CDS = ColumnDataSource(data=source_dict_subtype_mrna)
#
#     return subtype_protein_CDS, subtype_mRNA_CDS

    # The data dictionary will look like this:
    # 'x': [(geneA, basal), (geneA, her2)....]
    # 'y': [mean for geneA x basal, mean for geneA x her2...]
    # 'upper' [mean+SEM for geneA x basal, mean+SEM for geneA x her2...]
    # 'lower' [mean-SEM for geneA x basal, mean-SEM for geneA x her2...]
    # 'legend_group = [5 subtypes * 4]