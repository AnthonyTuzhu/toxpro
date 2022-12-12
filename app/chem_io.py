from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors
import pandas as pd

from sklearn.preprocessing import StandardScaler

def calc_descriptors_from_frame(df: pd.DataFrame, scale=True) -> pd.DataFrame:
    """ calculates rdkit descriptors from a smiles.txt file """

    df['ROMol'] = [Chem.MolFromInchi(inchi) for inchi in df.inchi]

    calc = MoleculeDescriptors.MolecularDescriptorCalculator([desc[0] for desc in Descriptors.descList])

    X = pd.DataFrame([list(calc.CalcDescriptors(mol)) for mol in df['ROMol']],
                     columns=list(calc.GetDescriptorNames()),
                     index=df.compound_id)
    X = X.loc[X.notnull().all(1), :]

    if scale:
        X = pd.DataFrame(StandardScaler().fit_transform(X), index=X.index, columns=X.columns)
    return X