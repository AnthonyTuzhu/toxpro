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
