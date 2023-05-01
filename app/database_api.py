from flask import Blueprint, request, jsonify, session
from flask_login import login_required, current_user
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase

from app.db_models import db
from app.db_models import User, Dataset, Chemical, QSARModel

import sys
import pandas as pd

print(sys.executable)

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

