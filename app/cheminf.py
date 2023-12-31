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

from sqlalchemy import exc

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as skl_PCA

import plotly
import plotly.express as px
import plotly.graph_objs as go

import json

import app.chem_io as chem_io
import app.machine_learning as ml

# from app.db import get_db
from app.db_models import User, db

import itertools

bp = Blueprint('cheminf', __name__, url_prefix='/cheminf')


@bp.route('/curator', methods=['GET', 'POST'])
@login_required
def curator():
    """
    displays the homepage

    """
    from app.curator.curator import Curator

    current_user.has_qsar = any(d.qsar_models for d in current_user.datasets)
    if request.method == 'GET':
        return render_template('cheminf/curator.html', user_datasets=list(current_user.datasets), user=current_user)

    dataset_selection = request.form['dataset-selection'].strip()
    dup_selection = request.form['duplicate-selection'].strip()
    create_or_replace = request.form['create-or-replace'].strip()
    if create_or_replace == 'replace':
        replace = True
    else:
        replace = False

    try:
        current_user.launch_task('curate_chems',
                                 f'Curating {dataset_selection} chemicals',
                                 current_user.id,
                                 dataset_selection,
                                 dup_selection,
                                 replace
                                 )
        db.session.commit()
    except exc.OperationalError as err:
        flash("Failed to submit curation job, please try again", 'error')

    return redirect(url_for('cheminf.curator'))


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
        .filter(Dataset.dataset_name == dataset_selection) \
        .filter(Dataset.user_id == current_user.id).statement
    df = pd.read_sql(query_statement, db.session.connection())

    desc_set = ['MolWt', 'TPSA', 'NumRotatableBonds', 'NumHDonors', 'NumHAcceptors', 'MolLogP']

    descriptors = chem_io.calc_descriptors_from_frame(df, scale=True, desc_set=desc_set)

    skl_PCA_fit = skl_PCA(n_components=3).fit(descriptors)
    pca = pd.DataFrame(skl_PCA_fit.transform(descriptors),
                       columns=['PCA1', 'PCA2', 'PCA3'])

    pca['CMP_ID'] = df.set_index('compound_id').loc[descriptors.index].index
    pca['Activity'] = df.set_index('compound_id').loc[descriptors.index].activity.values

    pca['Activity'] = pca['Activity'].astype('category')

    fig = px.scatter_3d(pca,
                        x='PCA1',
                        y='PCA2',
                        z='PCA3',
                        color='Activity',
                        hover_name='CMP_ID',
                        hover_data=['Activity', 'PCA1', 'PCA2'],
                        labels={
                            "PCA1": "PC1 ({:.2%})".format(skl_PCA_fit.explained_variance_ratio_[0]),
                            "PCA2": "PC2 ({:.2%})".format(skl_PCA_fit.explained_variance_ratio_[1]),
                            "PCA3": "PC3 ({:.2%})".format(skl_PCA_fit.explained_variance_ratio_[2]),
                        },
                        height=800,  # default height
                        color_discrete_map={
                            1: "rgba(255, 0, 0, 0.5)",
                            0: "rgba(0, 0, 255, 0.5)"
                        }
                        )

    fig.update_layout(template='plotly_white',
                      scene=dict(aspectratio=dict(x=1, y=1, z=1))
                      )
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
    # fig.update_zaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_traces(marker=dict(size=6,
                                  opacity=0.5,
                                  # line=dict(width=2, color='DarkSlateGrey')
                                  ),
                      )

    camera = dict(
        eye=dict(x=0., y=2.5, z=0),
    )

    fig.update_layout(scene_camera=camera, scene_dragmode='orbit')

    # fig.update_layout(shapes=[
    #     # unfilled rectange
    #     go.layout.Shape(
    #         type="rect",
    # xref ="paper",
    # yref ="paper",
    # x0 = 0,
    # y0 = -0.1,
    # x1 = 1.05,
    # y1 = 1,
    # line = {"width": 1, "color": "black"})])

    # Create graphJSON
    pca_plot = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    return render_template('cheminf/PCA.html', user_datasets=list(current_user.datasets), pca_plot=pca_plot)


@bp.route('/QSAR-build', methods=('GET', 'POST'))
@login_required
def QSAR_build():
    """ Registers a new user.  checkRecaptcha() must return True to register user.


    """
    current_user.has_qsar = any(d.qsar_models for d in current_user.datasets)
    if request.method == 'GET':
        return render_template('cheminf/QSAR-build.html', user_datasets=list(current_user.datasets), user=current_user)

    dataset_selection = request.form.getlist('dataset-selection')
    desc_selection = request.form.getlist('descriptor-selection')
    alg_selection = request.form.getlist('algorithm-selection')
    type_selection = request.form.getlist('type-selection')
    name_list = []
    for element in itertools.product(*[dataset_selection, desc_selection, alg_selection, type_selection]):
        name_list.append('&'.join(element))

    # name = f'{dataset_selection}-{desc_selection}-{alg_selection}-{type_selection}'
    for name in name_list:
        name_string = name.split("&")
        dataset = name_string[0]
        desc = name_string[1]
        alg = name_string[2]
        type1 = name_string[3]

        try:
            current_user.launch_task('build_qsar',
                                     f'Building QSAR Model on {name}',
                                     current_user.id,
                                     dataset,
                                     desc,
                                     alg,
                                     type1
                                     )


        except exc.OperationalError as err:
            flash("Failed to submit model, please try again", 'error')

        db.session.commit()
    return redirect(url_for('cheminf.QSAR_build'))


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
    model_name = request.form.getlist('model-selection')
    output_type = request.form['output-type'].strip()
    activity_col = request.form['smiles-column'].strip() or 'SMILES'

    error = None

    if not sdfile:
        error = "No SDFile was attached."

    if sdfile and not sdfile.filename.rsplit('.', 1)[1] in ['csv', 'sdf']:
        error = "The file is not an SDF"

    if sdfile:
        compound_filename = secure_filename(sdfile.filename)

        user_uploaded_file = os.path.join(current_app.instance_path, compound_filename)
        name = ntpath.basename(user_uploaded_file).split('.')[0]

        sdfile.save(user_uploaded_file)

        if sdfile and sdfile.filename.rsplit('.', 1)[1] in ['csv']:
            mols_df = pd.read_csv(user_uploaded_file)
            PandasTools.AddMoleculeColumnToFrame(mols_df, smilesCol=activity_col)
            mols_df = mols_df[mols_df.ROMol.notnull()]

            os.remove(user_uploaded_file)

            mols_df['compound_id'] = ['mol_{}'.format(i) for i in range(mols_df.shape[0])]
            mols_df['inchi'] = [Chem.MolToInchi(mol) for mol in mols_df.ROMol]
            mols_df_trim = mols_df

            if mols_df.empty:
                error = 'No compounds in SDFile'

            if not error:
                for model in model_name:
                    qsar_model = QSARModel.query.filter_by(name=model).first()
                    sklearn_model = pickle.loads(qsar_model.sklearn_model)

                    X_predict = chem_io.get_desc(mols_df_trim, qsar_model.descriptors)
                    mols_df_trim = mols_df_trim.set_index('compound_id').loc[X_predict.index]
                    mols_df_trim[f'{model}_Predictions'] = sklearn_model.predict(X_predict)
                    mols_df_trim.reset_index().drop(columns=['compound_id'])
                    mols_df_trim['compound_id'] = ['mol_{}'.format(i) for i in range(mols_df.shape[0])]

                if output_type == 'CSV':
                    import io
                    mem = io.BytesIO()
                    mem.write(mols_df_trim.drop(['ROMol', 'compound_id'], axis=1).to_csv().encode())
                    mem.seek(0)
                    return flask.send_file(
                        mem,
                        as_attachment=True,
                        download_name=f"{name}_predicted.csv",
                        mimetype="text/plain",
                    )

                if output_type == 'SDF':
                    download_file = os.path.join(current_app.instance_path + f"{compound_filename.split('.')[0]}.sdf")
                    PandasTools.WriteSDF(mols_df_trim, download_file,
                                         properties=mols_df_trim.drop('ROMol', axis=1).columns)
                    return flask.send_file(download_file,
                                           as_attachment=True,
                                           download_name=f"{name}_predicted.sdf")

        if sdfile and sdfile.filename.rsplit('.', 1)[1] in ['sdf']:
            mols_df = PandasTools.LoadSDF(user_uploaded_file)

            os.remove(user_uploaded_file)

            mols_df['compound_id'] = ['mol_{}'.format(i) for i in range(mols_df.shape[0])]
            mols_df['inchi'] = [Chem.MolToInchi(mol) for mol in mols_df.ROMol]
            mols_df_trim = mols_df

            if mols_df.empty:
                error = 'No compounds in SDFile'

            if not error:
                for model in model_name:
                    qsar_model = QSARModel.query.filter_by(name=model).first()
                    sklearn_model = pickle.loads(qsar_model.sklearn_model)

                    X_predict = chem_io.get_desc(mols_df_trim, qsar_model.descriptors)
                    mols_df_trim = mols_df_trim.set_index('compound_id').loc[X_predict.index]
                    mols_df_trim[f'{model}_Predictions'] = sklearn_model.predict(X_predict)
                    mols_df_trim.reset_index().drop(columns=['compound_id', 'ID'])
                    mols_df_trim['compound_id'] = ['mol_{}'.format(i) for i in range(mols_df.shape[0])]

                if output_type == 'CSV':
                    import io
                    mem = io.BytesIO()
                    mem.write(mols_df_trim.drop(['ROMol', 'compound_id', 'ID'], axis=1).to_csv().encode())
                    mem.seek(0)
                    return flask.send_file(
                        mem,
                        as_attachment=True,
                        download_name=f"{name}_predicted.csv",
                        mimetype="text/plain",
                    )

                if output_type == 'SDF':
                    PandasTools.WriteSDF(mols_df_trim, user_uploaded_file,
                                         properties=mols_df_trim.drop('ROMol', axis=1).columns)
                    return flask.send_file(user_uploaded_file,
                                           as_attachment=True,
                                           download_name=f"{name}_predicted.sdf")

    if error:
        flash(error, 'danger')

    return render_template('cheminf/QSAR-predict.html', user_qsar_models=user_qsar_models)


if __name__ == '__main__':
    from app import create_app

    app = create_app()
    # app.app_context().push()
    query_statement = db.session.query(Dataset).innerjoin(Dataset,
                                                          Dataset.id == Chemical.dataset_id) \
        .filter_by(dataset_name=1, user_id=1) \
        .first()
    print(query_statement)
