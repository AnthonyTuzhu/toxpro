from flask import Blueprint, request, jsonify, session
from flask_login import login_required, current_user
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase

from app.db_models import db
import app.config as config
from app.db_models import User, Dataset, Chemical, QSARModel

from app.master_db import get_database, make_query, get_master

import sys, os
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

@bp.route('/tox-bioprofile')
def get_bioprofile_data():
    endpoint_selection = request.args.get("endpointSelection", None)
    ep = TOXICITY_ENDPOINT_INFO.set_index('Endpoint').loc[endpoint_selection].to_dict()
    dataset = ep['Dataset']

    def confusion_matrix(df, activity_class, dataset_name):
        """ this function calculates the confusion matrix for an assay, toxicity pair """
        df[activity_class] = pd.to_numeric(df[activity_class], errors='coerce')
        df = df[df[activity_class].notnull()]

        tps = ((df[activity_class] == 1) & (df.Activity_Transformed == 1)).sum()
        fps = ((df[activity_class] == 0) & (df.Activity_Transformed == 1)).sum()
        tns = ((df[activity_class] == 0) & (df.Activity_Transformed == 0)).sum()
        fns = ((df[activity_class] == 1) & (df.Activity_Transformed == 0)).sum()

        return tps, fps, tns, fns

    bioprofile = pd.read_csv(os.path.join(config.Config.BIOPROFILE_DIR, f"{dataset}+{endpoint_selection}.csv"))
    med = (
        bioprofile[['Master-ID', endpoint_selection]]
        .drop_duplicates()
        [endpoint_selection]
        .median()
    )
    bioprofile['activity'] = bioprofile[endpoint_selection].copy()
    if dataset in ['LD50_curated']:
        bioprofile.loc[bioprofile[endpoint_selection] < med, 'activity'] = 1
        bioprofile.loc[bioprofile[endpoint_selection] >= med, 'activity'] = 0

    matrix = (
        bioprofile
        .groupby('AID')
        .apply(lambda x: confusion_matrix(x, endpoint_selection, dataset))
        .apply(pd.Series)
        .set_axis(['TP', 'FP', 'TN', 'FN'], axis=1)
        .reset_index()
        .sort_values('TP', ascending=False)
    )
    matrix['PPV'] = matrix.TP / (matrix.TP + matrix.FP)
    matrix['Sensitivity'] = matrix.TP / (matrix.TP + matrix.FN)
    bioprofile = pd.merge(matrix, bioprofile, on='AID', how='inner')

    PCA_DF = make_query('select [Master-ID], PCA1, PCA2, PCA3'
                        ' from chemical_space')

    CID_DF = make_query('select [Master-ID], [CID]'
                        ' from cid_lookup')

    pca = PCA_DF.merge(CID_DF, on='Master-ID', how='inner').join(
        bioprofile[['CID', 'activity']].drop_duplicates().set_index('CID'))
    pca['CIDs'] = pca.index
    pca = pd.merge(pca, bioprofile, on='Master-ID', how='inner')

    bio_info = pd.read_table(config.Config.BIOASSAYS)
    biodict = dict(zip(bio_info['AID'], bio_info['BioAssay Name']))
    pca['BioAssay Name'] = pca['AID'].map(biodict)

    table = pca.groupby(['AID', 'BioAssay Name'])['Activity_Transformed']\
        .value_counts()\
        .unstack(fill_value=0)\
        .rename(columns={-1.0: 'Inactive', 0.0: 'Inconclusive', 1.0: 'Active'})
    table['Active rate'] = table['Active'] / (table['Active'] + table['Inactive'])
    table = pd.DataFrame(table.to_records())
    top_assays = table.sort_values(by=['Active rate'], ascending=False).head(20).to_dict(orient='records')

    return top_assays