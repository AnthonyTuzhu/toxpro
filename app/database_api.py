from flask import Blueprint, request, jsonify, session
from flask_login import login_required, current_user
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase

from app.db_models import db
from app.db_models import User, Dataset, Chemical, QSARModel

from app.master_db import get_database, make_query, get_master

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
        'type': 'histogram',
        'marker': {
            'color': 'rgb(158,202,225)',
            'opacity': 0.6,
            'line': {
                'color': 'rgb(8,48,107)',
                'width': 1.5
            }
    }}

    return trace

@bp.route('/tox-pca-data')
def get_pca_data():
    """ send back two traces of data, with a specific color activity highlighted """

    endpoint_selection = request.args.get("endpointSelection", None)





    PCA_DF = make_query('select [Master-ID], PCA1, PCA2, PCA3'
                        ' from chemical_space')



    # this df contains the chemicals with the specific
    # endpoint we are looking for
    ep = TOXICITY_ENDPOINT_INFO.set_index('Endpoint').loc[endpoint_selection].to_dict()

    df = get_database(ep['Dataset'])
    df = df[['Master-ID', 'CleanedInChI', endpoint_selection]]
    df[endpoint_selection] = df[endpoint_selection].apply(pd.to_numeric, errors='coerce')
    df = df.dropna(subset=endpoint_selection)

    PCA_DF = PCA_DF.merge(df, how='left')

    PCA_DF.loc[PCA_DF[endpoint_selection].isnull(), 'color'] = 'rgba(20, 40, 186, 0.7)'
    PCA_DF.loc[PCA_DF[endpoint_selection].notnull(), 'color'] = 'rgb(201, 40, 0, 1)'

    # PCA_DF.loc[PCA_DF[endpoint_selection].isnull(), 'opacity'] = 0.2
    # PCA_DF.loc[PCA_DF[endpoint_selection].notnull(), 'opacity'] = 1

    PCA_DF.loc[PCA_DF[endpoint_selection].isnull(), 'size'] = 1
    PCA_DF.loc[PCA_DF[endpoint_selection].notnull(), 'size'] = 5

    trace = {
        'x': PCA_DF.PCA1.values.tolist(),
        'y': PCA_DF.PCA2.values.tolist(),
        'z': PCA_DF.PCA3.values.tolist(),
        'mode': 'markers',
        'marker': {
            'size': PCA_DF['size'].values.tolist(),
            'color': PCA_DF.color.values.tolist(),
            #'opacity': PCA_DF.opacity.values.tolist(),
            'line': {
                'width': 0
            },
        },
        'type': 'scatter3d'
    }
    return trace


@bp.route('/update-pca')
def update_pca():
    """ send back two traces of data, with a specific color activity highlighted """


    endpoint_selection = request.args.get("endpointSelection", None)

    PCA_DF = make_query('select [Master-ID], PCA1, PCA2, PCA3'
                        ' from chemical_space')

    # this df contains the chemicals with the specific
    # endpoint we are looking for
    ep = TOXICITY_ENDPOINT_INFO.set_index('Endpoint').loc[endpoint_selection].to_dict()

    df = get_database(ep['Dataset'])
    df = df[['Master-ID', 'CleanedInChI', endpoint_selection]]
    df[endpoint_selection] = df[endpoint_selection].apply(pd.to_numeric, errors='coerce')
    df = df.dropna(subset=endpoint_selection)

    PCA_DF = PCA_DF.merge(df, how='left')

    PCA_DF.loc[PCA_DF[endpoint_selection].isnull(), 'color'] = 'rgba(20, 40, 186, 0.7)'
    PCA_DF.loc[PCA_DF[endpoint_selection].notnull(), 'color'] = 'rgb(201, 40, 0, 1)'

    # PCA_DF.loc[PCA_DF[endpoint_selection].isnull(), 'opacity'] = 0.2
    # PCA_DF.loc[PCA_DF[endpoint_selection].notnull(), 'opacity'] = 1

    PCA_DF.loc[PCA_DF[endpoint_selection].isnull(), 'size'] = 1
    PCA_DF.loc[PCA_DF[endpoint_selection].notnull(), 'size'] = 5

    update = {
        'marker.color':  [PCA_DF.color.values.tolist()],
        'marker.size': [PCA_DF['size'].values.tolist()]
    }
    return update

