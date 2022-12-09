# this module is outlined
# here: https://flask.palletsprojects.com/en/2.0.x/tutorial/views/


from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, current_app
)
from flask_login import current_user
from app.db_models import User, Dataset, Chemical
import pandas as pd

from rdkit import Chem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as skl_PCA

import plotly
import plotly.express as px

import json

#from app.db import get_db
from app.db_models import User, db

bp = Blueprint('cheminf', __name__, url_prefix='/cheminf')

@bp.route('/PCA', methods=('GET', 'POST'))
def PCA():
    """ Registers a new user.  checkRecaptcha() must return True to register user.


    """
    if request.method == 'GET':
        return render_template('cheminf/PCA.html', user_datasets=current_user.datasets)

    dataset_selection = request.form['dataset-selection'].strip()

    query_statement = db.session.query(Chemical).join(Dataset,
                                                       Dataset.id == Chemical.dataset_id) \
                        .filter(Dataset.dataset_name==dataset_selection) \
                        .filter(Dataset.user_id==current_user.id).statement
    df = pd.read_sql(query_statement, db.session.bind)
    df['ROMol'] = [Chem.MolFromInchi(inchi) for inchi in df.inchi]

    descriptors = calc_descriptors_from_frame(df)

    print(df.columns)

    pca = pd.DataFrame(skl_PCA(n_components=2).fit_transform(descriptors),
                       columns=['PCA1', 'PCA2'])

    pca['CMP_ID'] = df.set_index('compound_id').loc[descriptors.index].index
    pca['Activity'] = df.set_index('compound_id').loc[descriptors.index].activity.values

    print(pca)
    fig = px.scatter(pca,
                     x='PCA1',
                     y='PCA2',
                     color='Activity'
                     )

    # Create graphJSON
    pca_plot = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    return render_template('cheminf/PCA.html', user_datasets=current_user.datasets, pca_plot=pca_plot)


@bp.route('/QSAR-build', methods=('GET', 'POST'))
def QSAR_build():
    """ Registers a new user.  checkRecaptcha() must return True to register user.


    """
    if request.method == 'GET':
        return render_template('cheminf/QSAR-build.html')

    return render_template('cheminf/QSAR-build.html')


@bp.route('/QSAR-predict', methods=('GET', 'POST'))
def QSAR_predict():
    """ Registers a new user.  checkRecaptcha() must return True to register user.


    """
    if request.method == 'GET':
        return render_template('cheminf/QSAR-build.html')

    return render_template('cheminf/QSAR-build.html')



def calc_descriptors_from_frame(df: pd.DataFrame) -> pd.DataFrame:
    """ calculates rdkit descriptors from a smiles.txt file """


    calc = MoleculeDescriptors.MolecularDescriptorCalculator([desc[0] for desc in Descriptors.descList])

    X = pd.DataFrame([list(calc.CalcDescriptors(mol)) for mol in df['ROMol']],
                     columns=list(calc.GetDescriptorNames()),
                     index=df.compound_id)
    X = X.loc[X.notnull().all(1), :]

    X = pd.DataFrame(StandardScaler().fit_transform(X), index=X.index, columns=X.columns)
    return X
if __name__ == '__main__':
    from app import create_app

    app = create_app()
    #app.app_context().push()
    query_statement = db.session.query(Dataset).innerjoin(Dataset,
                                                       Dataset.id == Chemical.dataset_id) \
                        .filter_by(dataset_name=1, user_id=1) \
                        .first()
    print(query_statement)