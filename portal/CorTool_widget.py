import sys
import pandas as pd


def protein_import_data():
    """read in lower triangle correlation matrix generated from protein data"""
    jo_pro_cor_df = pd.read_csv("Data/Johansson/jo_pro_tricormat.csv", index_col="Gene")
    kr_pro_cor_df = pd.read_csv("Data/Krug/kr_pro_tricormat.csv", index_col="Gene")
    me_pro_cor_df = pd.read_csv("Data/Mertins/me_pro_tricormat.csv", index_col="Gene")

    return jo_pro_cor_df, kr_pro_cor_df, me_pro_cor_df


def mrna_import_data():
    """read in lower triangle correlation matrix generated from mrna data"""
    jo_mrna_cor_df = pd.read_csv("Data/Johansson/jo_mrna_tricormat.csv", index_col="Gene")
    kr_mrna_cor_df = pd.read_csv("Data/Krug/kr_mrna_tricormat.csv", index_col="Gene")
    me_mrna_cor_df = pd.read_csv("Data/Mertins/me_mrna_tricormat.csv", index_col="Gene")

    return jo_mrna_cor_df, kr_mrna_cor_df, me_mrna_cor_df


def df_to_dict_converter(df, gene) -> dict:
    """convert the df into dict in form like this {'index': top 20 related gene names, 'data': correlation values}"""
    df = df[gene].sort_values(ascending=True).head(20)
    df = pd.DataFrame(df, columns=None)
    dictionary = df.to_dict('split')
    del dictionary['columns']
    return dictionary


def button_callback():
    """call back function to stop the server"""
    sys.exit()  # Stop the server
