# this module is outlined
# here: https://flask.palletsprojects.com/en/2.0.x/tutorial/views/
import flask
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, current_app
)
from flask_login import current_user, login_required
from app.db_models import User, Dataset, Chemical, QSARModel
import pandas as pd
import ntpath, os

from werkzeug.utils import secure_filename
from rdkit import Chem
from rdkit.Chem import PandasTools
import pickle


from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as skl_PCA

import plotly
import plotly.express as px

import json

import app.chem_io as chem_io
import app.machine_learning as ml

#from app.db import get_db
from app.db_models import User, db

bp = Blueprint('cheminf', __name__, url_prefix='/cheminf')

@bp.route('/PCA', methods=('GET', 'POST'))
@login_required
def PCA():
    """ Registers a new user.  checkRecaptcha() must return True to register user.


    """
    if request.method == 'GET':
        return render_template('cheminf/PCA.html', user_datasets=list(current_user.datasets))

    dataset_selection = request.form['dataset-selection'].strip()

    query_statement = db.session.query(Chemical).join(Dataset,
                                                       Dataset.id == Chemical.dataset_id) \
                        .filter(Dataset.dataset_name==dataset_selection) \
                        .filter(Dataset.user_id==current_user.id).statement
    df = pd.read_sql(query_statement, db.session.bind)


    descriptors = chem_io.calc_descriptors_from_frame(df)

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

    return render_template('cheminf/PCA.html', user_datasets=list(current_user.datasets), pca_plot=pca_plot)


@bp.route('/QSAR-build', methods=('GET', 'POST'))
@login_required
def QSAR_build():
    """ Registers a new user.  checkRecaptcha() must return True to register user.


    """

    if request.method == 'GET':
        return render_template('cheminf/QSAR-build.html', user_datasets=list(current_user.datasets))


    dataset_selection = request.form['dataset-selection'].strip()
    desc_selection = request.form['descriptor-selection'].strip()
    alg_selection = request.form['algorithm-selection'].strip()


    query_statement = db.session.query(Chemical).join(Dataset,
                                                       Dataset.id == Chemical.dataset_id) \
                        .filter(Dataset.dataset_name==dataset_selection) \
                        .filter(Dataset.user_id==current_user.id).statement
    df = pd.read_sql(query_statement, db.session.bind)


    dataset = Dataset.query.filter_by(dataset_name=dataset_selection).first()
    # erase exists models
    qsar_model = QSARModel.query.filter_by(user_id=current_user.id,
                                       algorithm=alg_selection,
                                       descriptors=desc_selection,
                                       dataset_id=dataset.id,).first()

    if qsar_model:
        db.session.delete(qsar_model)


    # create descriptors
    X = chem_io.get_desc(df, desc_selection)

    y = df['activity']
    y.index = df['compound_id']
    y = y.loc[X.index]

    if desc_selection == 'RDKit':
        scale = True
    else:
        scale = False

    model, cv_preds, train_stats = ml.build_qsar_model(X,
                                                       y,
                                                       alg_selection,
                                                       scale=scale)


    name = f'{dataset_selection}-{desc_selection}-{alg_selection}'



    qsar_model = QSARModel(user_id=current_user.id,
                           name=name,
                           algorithm=alg_selection,
                           descriptors=desc_selection,
                           dataset_id=dataset.id,
                           sklearn_model=model)

    db.session.add(qsar_model)
    db.session.commit()

    flash(f'Finished QSAR model using {desc_selection} '
          f'descriptors and {alg_selection} on {dataset_selection}.  Model saved as {qsar_model.name}', 'info')

    qsar_results=train_stats
    print(qsar_results)
    return render_template('cheminf/QSAR-build.html',
                           user_datasets=list(current_user.datasets),
                           qsar_results=qsar_results)


@bp.route('/QSAR-predict', methods=('GET', 'POST'))
@login_required
def QSAR_predict():
    """ Registers a new user.  checkRecaptcha() must return True to register user.


    """


    user_qsar_models = list(QSARModel.query.with_entities(QSARModel.user_id,
                                                          QSARModel.name
                                                          ).filter_by(user_id=current_user.id).all())

    if request.method == 'GET':
        return render_template('cheminf/QSAR-predict.html', user_qsar_models=user_qsar_models)


    sdfile = request.files['predict-file']
    model_name = request.form['model-selection'].strip()


    error = None

    if not sdfile:
        error = "No SDFile was attached."

    if sdfile and not sdfile.filename.rsplit('.', 1)[1] in ['sdf']:
        error = "The file is not an SDF"

    if sdfile:
        compound_filename = secure_filename(sdfile.filename)

        user_uploaded_file = os.path.join(current_app.instance_path, compound_filename)
        name = ntpath.basename(user_uploaded_file).split('.')[0]


        sdfile.save(user_uploaded_file)

        mols_df = PandasTools.LoadSDF(user_uploaded_file)
        os.remove(user_uploaded_file)

        mols_df['compound_id'] = ['mol_{}'.format(i) for i in range(mols_df.shape[0])]
        mols_df['inchi'] = [Chem.MolToInchi(mol) for mol in mols_df.ROMol]

        if mols_df.empty:
            error = 'No compounds in SDFile'

        if not error:

            qsar_model = QSARModel.query.filter_by(name=model_name).first()
            sklearn_model = pickle.loads(qsar_model.sklearn_model)

            X_predict = chem_io.get_desc(mols_df, qsar_model.descriptors)

            mols_df_trim = mols_df.set_index('compound_id').loc[X_predict.index]

            mols_df_trim[f'{model_name}_Predictions'] = sklearn_model.predict(X_predict)

            PandasTools.WriteSDF(mols_df_trim, user_uploaded_file,
                                 properties=mols_df_trim.drop('ROMol', axis=1).columns)
            return flask.send_file(user_uploaded_file,
                                   attachment_filename=compound_filename.replace('.sdf', '_predicted.sdf'))


    if error:
        flash(error, 'danger')

    return render_template('cheminf/QSAR-predict.html', user_qsar_models=user_qsar_models)


if __name__ == '__main__':
    from app import create_app

    app = create_app()
    #app.app_context().push()
    query_statement = db.session.query(Dataset).innerjoin(Dataset,
                                                       Dataset.id == Chemical.dataset_id) \
                        .filter_by(dataset_name=1, user_id=1) \
                        .first()
    print(query_statement)