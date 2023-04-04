import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import PandasTools
import config, os
import glob, ntpath

DATA_DIR = config.Config.DATA_DIR
# get names of all the databases in
# the data directory folder
ALL_DATABASES = [ntpath.basename(f)[:-4] for f in glob.glob(os.path.join(DATA_DIR, '*.sdf'))]

def load_db(db_name: str, mol_col=False) -> pd.DataFrame:
    """ load a specific dataset as a pandas dataframe

    db_name: string name of the database
    mol_col: bool, whether to load the rdkit mol columns
    """
    master_sdf = os.path.join(DATA_DIR, '{}.sdf'.format(db_name))

    df = PandasTools.LoadSDF(master_sdf).replace('nan', np.nan)
    if not mol_col and 'ROMol' in df.columns:
        df = df.drop('ROMol', axis=1)
    return df

def load_master(mol_col=False) -> pd.DataFrame:
    return load_db('Merged', mol_col=mol_col)

def print_overview():
    df = load_master()
    print("There are {} compounds in the master database.".format(df.shape[0]))
    print("Here is an overview of the databases:")
    for db_name in ALL_DATABASES:
        if db_name == 'Merged':
            continue
        df = load_db(db_name)
        print("The database {} has {} compounds.".format(db_name, df.shape[0]))


if __name__ == '__main__':
    print_overview()


