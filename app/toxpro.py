from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, jsonify, send_from_directory, current_app,
    send_file, abort, Response
)
from werkzeug.utils import secure_filename

from rdkit import Chem
from rdkit.Chem import PandasTools

import plotly
import plotly.express as px
import json, os, ntpath

from app.db_models import User, Dataset, Chemical, db, Permission
import app.master_db as master_db
from flask_login import current_user, login_required, login_manager
from sqlalchemy import exc
import pandas as pd
from plotly.graph_objs import *
from functools import wraps

import app.pubchem as pc

bp = Blueprint('toxpro', __name__)

TOXICITY_ENDPOINT_INFO = pd.read_csv('data/toxicity-endpoint-info.csv', index_col=0)


def permission_required(permission):
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            if not current_user.can(permission):
                abort(403)
                return f(*args, **kwargs)

        return decorated_function()

    return decorator


def advanced_required(f):
    return permission_required(Permission.ADVANCED)(f)


# this is necessary for declaring
# variables that are available across
# all templates: https://stackoverflow.com/questions/26498689/flask-jinja-pass-data-to-a-base-template-all-templates
# @bp.context_processor
# def inject_nms():
#     #

@bp.route('/', methods=['GET'])
def index():
    """
    displays the homepage

    """
    return render_template('toxpro/home.html')


@bp.route('/about', methods=['GET'])
def about():
    """
    displays the about page

    """
    return render_template('toxpro/about.html')


@bp.route('/contact', methods=['GET'])
def contact():
    """
    displays the contact homepage

    """
    return render_template('toxpro/contact.html')


@bp.route('/datasets', methods=['GET'])
@login_required
def datasets():
    """
    displays the homepage

    """
    return render_template('toxpro/datasets.html', user_datasets=list(current_user.datasets))


@bp.route('/AdvancedFunctions', methods=['GET'])
@login_required
def AdvancedFunctions():
    """
    displays the homepage

    """
    return render_template('toxpro/advancedfunctions.html')


@bp.route('/upload_dataset', methods=['POST'])
@login_required
def upload_dataset():
    """
    uploads a dataset

    For the javascript datatable on the client-side, there is a nice
    tutorial overview here: https://blog.miguelgrinberg.com/post/beautiful-interactive-tables-for-your-flask-templates/page/0


    his basic pagination tutorial is here: https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-ix-pagination

    Both are necessary for the
    """
    sdfile = request.files['compound_file']
    activity_col = request.form['act-id-property'].strip() or 'Activity'
    compound_id_col = request.form['cmp-id-property'].strip() or 'CMP_ID'
    smiles_col = request.form['smiles-id-property'].strip() or 'SMILES'
    dataset_type = request.form['dataset-type'].strip()

    error = None

    if not sdfile:
        error = "No SDFile was attached."

    if sdfile and not sdfile.filename.rsplit('.', 1)[1] in ['csv','sdf']:
        error = "The file is not an SDF"

    if sdfile:
        compound_filename = secure_filename(sdfile.filename)

        user_uploaded_file = os.path.join(current_app.instance_path, compound_filename)
        name = ntpath.basename(user_uploaded_file).split('.')[0]

        # create the dataset
        try:
            dataset = Dataset(user_id=current_user.id, dataset_name=name, type=dataset_type)
            db.session.add(dataset)
            db.session.commit()
        except exc.IntegrityError:
            db.session.rollback()
            error = f'Sorry, there is already a dataset named {name}'
            flash(error, 'danger')
            return redirect(url_for('toxpro.datasets'))

        # I think we have to save this in order to use it, not sure if we can read it otherwise
        sdfile.save(user_uploaded_file)

        if sdfile and sdfile.filename.rsplit('.', 1)[1] in ['csv']:
            mols_df = pd.read_csv(user_uploaded_file)
            PandasTools.AddMoleculeColumnToFrame(mols_df, smilesCol=smiles_col)
            mols_df = mols_df[mols_df.ROMol.notnull()]

        if sdfile and sdfile.filename.rsplit('.', 1)[1] in ['sdf']:
            mols_df = PandasTools.LoadSDF(user_uploaded_file)

        os.remove(user_uploaded_file)

        if mols_df.empty:
            error = 'No compounds in the file'
        if activity_col not in mols_df.columns:
            error = f'Activity {activity_col} not in the file.'
        if compound_id_col not in mols_df.columns:
            error = f'Compound ID {compound_id_col} not in the file.'

        if error == None and dataset_type == 'Binary':

            # coerce activity column to be
            # integer
            mols_df[activity_col] = pd.to_numeric(mols_df[activity_col], errors='coerce')
            mols_df = mols_df[mols_df[activity_col].notnull()]
            mols_df = mols_df[(mols_df[activity_col] == 0) | (mols_df[activity_col] == 1)]
            mols_df[activity_col] = mols_df[activity_col].astype(int)

            mols_df = mols_df.sort_values(activity_col, ascending=False)
            mols_df = mols_df.drop_duplicates(compound_id_col, keep='first')

            for i, row in mols_df.iterrows():
                mol = row.ROMol
                activity = row[activity_col]
                cmp_id = row[compound_id_col]
                inchi = Chem.MolToInchi(mol)

                if cmp_id and inchi and (activity in [1, 0, '1', '0', '1.0', '0.0']):
                    chem = Chemical(inchi=inchi, dataset_id=dataset.id, activity=activity, compound_id=cmp_id)
                    dataset.chemicals.append(chem)

            db.session.add(dataset)
            db.session.commit()

            num_chemicals = len(list(dataset.chemicals))
            flash(f'Uploaded {name} as a new dataset; Added {num_chemicals} chemicals', 'success')
            return redirect(url_for('toxpro.datasets'))

        elif error == None and dataset_type == 'Continuous':
            mols_df[activity_col] = pd.to_numeric(mols_df[activity_col], errors='coerce')
            mols_df = mols_df[mols_df[activity_col].notnull()]
            mols_df[activity_col] = mols_df[activity_col].astype(int)

            mols_df = mols_df.sort_values(activity_col, ascending=False)
            mols_df = mols_df.drop_duplicates(compound_id_col, keep='first')

            for i, row in mols_df.iterrows():
                mol = row.ROMol
                activity = row[activity_col]
                cmp_id = row[compound_id_col]
                inchi = Chem.MolToInchi(mol)

                if cmp_id and inchi and activity:
                    chem = Chemical(inchi=inchi, dataset_id=dataset.id, activity=activity, compound_id=cmp_id)
                    dataset.chemicals.append(chem)

            db.session.add(dataset)
            db.session.commit()

            num_chemicals = len(list(dataset.chemicals))
            flash(f'Uploaded {name} as a new dataset; Added {num_chemicals} chemicals', 'success')
            return redirect(url_for('toxpro.datasets'))

        else:
            db.session.delete(dataset)
            db.session.commit()

    flash(error, 'danger')

    return redirect(url_for('toxpro.datasets'))


@bp.route('/remove_dataset', methods=['POST'])
@login_required
def remove_dataset():
    """
    uploads a dataset

    """

    dataset_selection = request.form['dataset-selection'].strip()
    do_what_with_dataset = request.form['action']

    # there are two buttons one to download and one to remove.
    if do_what_with_dataset == 'Download dataset as CSV file':
        query_statement = db.session.query(Chemical).join(Dataset,
                                                          Dataset.id == Chemical.dataset_id) \
            .filter(Dataset.dataset_name == dataset_selection) \
            .filter(Dataset.user_id == current_user.id).statement

        df = pd.read_sql(query_statement, db.session.connection())

        import io
        mem = io.BytesIO()
        mem.write(df.to_csv().encode())
        mem.seek(0)
        return send_file(
            mem,
            as_attachment=True,
            download_name=f"{dataset_selection}.csv",
            mimetype="text/plain",
        )

    dataset = Dataset.query.filter_by(dataset_name=dataset_selection, user_id=current_user.id).first()
    db.session.delete(dataset)
    db.session.commit()
    message = f"Removed {dataset_selection}"
    flash(message, 'danger')

    return redirect(url_for('toxpro.datasets'))


@bp.route('/import_pubchem', methods=['POST'])
@login_required
def import_pubchem():
    """
    import a pubchem data set
    """

    pubchem_aid_string = request.form.get('pubchem_aid', None)

    try:
        pubchem_aid = int(pubchem_aid_string)
    except ValueError:
        pubchem_aid = None

    if pubchem_aid == None:
        flash(f'"{pubchem_aid_string}" is not a valid PubChem AID', 'danger')
        return redirect(url_for('toxpro.datasets'))

    name, reason = pc.get_assay_name(pubchem_aid)
    if reason is not None:
        flash(f'Could not find AID {pubchem_aid}', 'danger')
        return redirect(url_for('toxpro.datasets'))

    current_user.launch_task('add_pubchem_data',
                             f'Importing structure-activity data for AID {pubchem_aid}: {name}',
                             pubchem_aid,
                             current_user.id
                             )
    db.session.commit()
    flash(f'Successfully submitted job for gathering structure-activity data for AID {pubchem_aid}: {name}', 'success')
    return redirect(url_for('toxpro.datasets'))
    # if not sdfile:
    #     error = "No SDFile was attached."
    #
    # if sdfile and not sdfile.filename.rsplit('.', 1)[1] in ['sdf']:
    #     error = "The file is not an SDF"
    #
    # if sdfile:
    #     compound_filename = secure_filename(sdfile.filename)
    #
    #     user_uploaded_file = os.path.join(current_app.instance_path, compound_filename)
    #     name = ntpath.basename(user_uploaded_file).split('.')[0]
    #
    #     # create the dataset
    #     try:
    #         dataset = Dataset(user_id=current_user.id, dataset_name=name)
    #         db.session.add(dataset)
    #         db.session.commit()
    #     except exc.IntegrityError:
    #         db.session.rollback()
    #         error = f'Sorry, there is already a dataset named {name}'
    #         flash(error, 'danger')
    #         return redirect(url_for('toxpro.datasets'))
    #
    #     # I think we have to save this in order to use it, not sure if we car read it otherwise
    #     sdfile.save(user_uploaded_file)
    #
    #     mols_df = PandasTools.LoadSDF(user_uploaded_file)
    #     os.remove(user_uploaded_file)
    #
    #     if mols_df.empty:
    #         error = 'No compounds in SDFile'
    #     if activity_col not in mols_df.columns:
    #         error = f'Activity {activity_col} not in SDFile.'
    #     if compound_id_col not in mols_df.columns:
    #         error = f'Compound ID {compound_id_col} not in SDFile.'
    #
    #     if error == None:
    #
    #         # coerce activity column to be
    #         # integer
    #         mols_df[activity_col] = pd.to_numeric(mols_df[activity_col], errors='coerce')
    #         mols_df = mols_df[mols_df[activity_col].notnull()]
    #         mols_df = mols_df[(mols_df[activity_col] == 0) | (mols_df[activity_col] == 1)]
    #         mols_df[activity_col] = mols_df[activity_col].astype(int)
    #
    #         mols_df = mols_df.sort_values(activity_col, ascending=False)
    #         mols_df = mols_df.drop_duplicates(compound_id_col, keep='first')
    #
    #         for i, row in mols_df.iterrows():
    #             mol = row.ROMol
    #             activity = row[activity_col]
    #             cmp_id = row[compound_id_col]
    #             inchi = Chem.MolToInchi(mol)
    #
    #             if cmp_id and inchi and (activity in [1, 0, '1', '0', '1.0', '0.0']):
    #                 chem = Chemical(inchi=inchi, dataset_id=dataset.id, activity=activity, compound_id=cmp_id)
    #                 dataset.chemicals.append(chem)
    #
    #         db.session.add(dataset)
    #         db.session.commit()
    #
    #         num_chemicals = len(list(dataset.chemicals))
    #         flash(f'Uploaded {name} as a new dataset; Added {num_chemicals} chemicals', 'success')
    #         return redirect(url_for('toxpro.datasets'))
    #     else:
    #         db.session.delete(dataset)
    #         db.session.commit()
    #
    # flash(error, 'danger')

    return redirect(url_for('toxpro.datasets'))


@bp.route('/assayProfile', methods=['GET'])
@login_required
def assayProfile():
    """
    displays the homepage

    """

    return render_template('toxpro/profile.html')


@bp.route('/toxdata', methods=['GET'])
@login_required
def toxdata():
    return render_template('toxpro/toxdata.html',
                           current_dbs=master_db.CURRENT_DATABASES,
                           endpoints=TOXICITY_ENDPOINT_INFO.to_dict('records'))


@bp.route('/download_database', methods=['POST'])
@login_required
def download_database():
    database_selection = request.form['database-selection'].strip()
    db = master_db.get_database(database_selection)

    import io
    mem = io.BytesIO()
    mem.write(db.to_csv().encode())
    mem.seek(0)
    return send_file(
        mem,
        as_attachment=True,
        download_name=f"{database_selection}.csv",
        mimetype="text/plain",
    )
