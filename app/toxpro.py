from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, jsonify, send_from_directory, current_app
)
from werkzeug.utils import secure_filename

from rdkit import Chem
from rdkit.Chem import PandasTools

import plotly
import plotly.express as px
import json, os, ntpath

from app.db_models import User, Dataset, Chemical, db
from flask_login import current_user, login_required
from sqlalchemy import exc

bp = Blueprint('toxpro', __name__)


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


@bp.route('/upload_dataset', methods=['POST'])
@login_required
def upload_dataset():
    """
    uploads a dataset

    """
    sdfile = request.files['compound_file']
    activity_col = request.form['act-id-property'].strip() or 'Activity'
    compound_id_col = request.form['cmp-id-property'].strip() or 'CMP_ID'

    error = None

    if not sdfile:
        error = "No SDFile was attached."

    if sdfile and not sdfile.filename.rsplit('.', 1)[1] in ['sdf']:
        error = "The file is not an SDF"

    if sdfile:
        compound_filename = secure_filename(sdfile.filename)

        user_uploaded_file = os.path.join(current_app.instance_path, compound_filename)
        name = ntpath.basename(user_uploaded_file).split('.')[0]

        # create the dataset
        try:
            dataset = Dataset(user_id=current_user.id, dataset_name=name)
            db.session.add(dataset)
            #db.session.commit()
        except exc.IntegrityError:
            db.session.rollback()
            error = f'Sorry, there is already a dataset named {name}'
            flash(error, 'danger')
            return redirect(url_for('toxpro.datasets'))

        # I think we have to save this in order to use it, not sure if we car read it otherwise
        sdfile.save(user_uploaded_file)

        mols_df = PandasTools.LoadSDF(user_uploaded_file)
        os.remove(user_uploaded_file)

        if mols_df.empty:
            error = 'No compounds in SDFile'
        if activity_col not in mols_df.columns:
            error = f'Activity {activity_col} not in SDFile.'
        if compound_id_col not in mols_df.columns:
            error = f'Compound ID {compound_id_col} not in SDFile.'

        if error == None:
            for i, row in mols_df.iterrows():
                mol = row.ROMol
                activity = row[activity_col]
                cmp_id = row[compound_id_col]
                inchi = Chem.MolToInchi(mol)

                if activity and cmp_id and inchi and (activity in [1, 0, '1', '0', '1.0', '0.0']):
                    chem = Chemical(inchi=inchi, dataset_id=dataset.id, activity=activity, compound_id=cmp_id)
                    dataset.chemicals.append(chem)

            db.session.add(dataset)
            db.session.commit()

            num_chemicals =len(list(dataset.chemicals))
            flash(f'Uploaded {name} as a new dataset; Added {num_chemicals} chemicals', 'success')
            return redirect(url_for('toxpro.datasets'))


    flash(error, 'danger')

    return redirect(url_for('toxpro.datasets'))


@bp.route('/remove_dataset', methods=['POST'])
@login_required
def remove_dataset():
    """
    uploads a dataset

    """

    dataset_selection = request.form['dataset-selection'].strip()

    dataset = Dataset.query.filter_by(dataset_name=dataset_selection).first()
    print(dataset)
    db.session.delete(dataset)
    db.session.commit()
    message = f"Removed {dataset_selection}"
    flash(message, 'danger')

    return redirect(url_for('toxpro.datasets'))


@bp.route('/curator', methods=['GET'])
@login_required
def curator():
    """
    displays the homepage

    """

    return render_template('toxpro/curator.html')


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
    """
    displays the homepage

    """

    return render_template('toxpro/toxdata.html')