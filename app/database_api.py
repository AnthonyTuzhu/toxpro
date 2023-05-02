from flask import Blueprint, request, jsonify, session
from flask_login import login_required, current_user
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase

from app.db_models import db
from app.db_models import User, Dataset, Chemical, QSARModel

from app.master_db import get_database

import sys
import pandas as pd

TOXICITY_ENDPOINT_INFO = pd.read_csv('data/toxicity-endpoint-info.csv', index_col=0)

bp = Blueprint('api', 'api', url_prefix='/api')

@bp.route('/datasets/<dataset_selection>')
def get_dataset_overview(dataset_selection):


    query_statement = db.session.query(Chemical).join(Dataset,
                                                       Dataset.id == Chemical.dataset_id) \
                                                        .filter(Dataset.dataset_name==dataset_selection) \
                                                        .filter(Dataset.user_id==current_user.id).statement
    df = pd.read_sql(query_statement, db.session.connection())

    actives = df['activity'].sum()
    inactives = df.shape[0] - actives


    data = {
        'name': dataset_selection,
        'actives': int(actives),
        'inactives': int(inactives)
    }
    return jsonify(results=data)

@bp.route('/dataset-data')
def get_dataset_data():

    dataset_selection = request.args.get("datasetSelection")
    search = request.args.get('search[value]')


    query = Dataset.query.filter(Dataset.dataset_name==dataset_selection) \
                                   .filter(Dataset.user_id==current_user.id) \
                                   .one().get_chemicals()


    total_chemicals = query.count()

    if search:
        query = query.filter(db.or_(
            Chemical.compound_id.like(f'%{search}%'),
            Chemical.activity.like(f'%{search}%')
        ))
        print(search)
    total_filtered = query.count()


    # pagination
    start = request.args.get('start', type=int)
    length = request.args.get('length', type=int)
    query = query.offset(start).limit(length)

    return {
        "data": [chemical.to_dict(structure_as_svg=True) for chemical in query.all()],
        'recordsFiltered': total_filtered,
        'recordsTotal': total_chemicals,
        'draw': request.args.get('draw', type=int)
    }

@bp.route('/tox-ep')
def get_toxicity_endpoint():
    """ send back a dataset endpoint """
    endpoint_selection = request.args.get("endpointSelection")
    ep = TOXICITY_ENDPOINT_INFO.set_index('Endpoint').loc[endpoint_selection].to_dict()

    df = get_database(ep['Dataset'])
    df = df[['Master-ID', 'CleanedInChI', endpoint_selection]]
    df[endpoint_selection] = df[endpoint_selection].apply(pd.to_numeric, errors='coerce')
    df = df.dropna(subset=endpoint_selection)
    trace = {
        'x': df[endpoint_selection].values.tolist(),
        'type': 'histogram'
    }
    return trace

