import sqlite3 as sql
import pandas as pd
import app.config as config
from rdkit.Chem import PandasTools

def make_connection():
    return sql.connect(config.Config.MASTER_DB_FILE)

def get_master() -> pd.DataFrame:
    con = make_connection()
    return pd.read_sql('select * from Master_database', con=con)

def make_query(query: str) -> pd.DataFrame:
    con = make_connection()
    return pd.read_sql(query, con=con)


def get_database(db_name: str) -> pd.DataFrame:
    """ gets all compounds in the database associated with a give table,
     the db_name should be the chemical table name """

    q = "select ml.[Master-ID] as [Master-ID], " \
        "ml.CleanedInChI as CleanedInChI, " \
        "cl.CID as CID_COL, " \
        "db.* " \
        "from master_lookup ml " \
        "inner join cid_lookup cl on ml.[Master-ID] = cl.[Master-ID] " \
        "inner join {} db on db.[Dataset-ID] = ml.[Dataset-ID] " \
        "WHERE db == '{}' ".format(db_name, db_name)

    return make_query(q)

TABLE_ACTIVITES = {
    'Hepatotoxicity_curated': ['PC_HT_class', 'H_HT_class'],
    'DART_curated': ['Oral-dev', 'Oral-mat'],
    'LD50_curated': ['LD50_mgkg'],
    'Estrogen_curated': ['ER-alpha log10 of 100*RBA (relative binding affinity vs estrogen)','Agonist Class','Agonist Potency','Antagonist Class','Antagonist Potency','Binding Class','Binding Potency','Uterotrophic Class']
}

CURRENT_DATABASES = [record for record in make_query('select distinct db from master_lookup').db.values
                     if '_curated' in record]

if __name__ == '__main__':
    hep = get_database('Hepatotoxicity_curated')
    print(hep)


