import sqlite3 as sql
import pandas as pd
import app.config as config

TABLE_ACTIVITES = {
    'Hepatotoxicity_curated': ['PC_HT_class', 'H_HT_class'],
    'DART_curated': ['Oral-dev', 'Oral-mat'],
    'LD50_curated': ['LD50_mgkg'],
    'Estrogen_curated': ['ER-alpha log10 of 100*RBA (relative binding affinity vs estrogen)','Agonist Class','Agonist Potency','Antagonist Class','Antagonist Potency','Binding Class','Binding Potency','Uterotrophic Class']
}


def make_connection():
    return sql.connect(config.Config.MASTER_DB_FILE)

def get_master() -> pd.DataFrame:
    con = make_connection()
    return pd.read_sql('select * from Master_database', con=con)

def make_query(query: str) -> pd.DataFrame:
    con = make_connection()
    return pd.read_sql(query, con=con)

