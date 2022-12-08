# this module is outlined
# here: https://flask.palletsprojects.com/en/2.0.x/tutorial/views/


from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for
)

#from app.db import get_db
from app.db_models import User, db

bp = Blueprint('cheminf', __name__, url_prefix='/cheminf')

@bp.route('/upload', methods=('GET', 'POST'))
def register():
    """ Registers a new user.  checkRecaptcha() must return True to register user.


    """
    if request.method == 'GET':
        return render_template('cheminf/upload.html')

    return redirect(url_for('cheminf.upload'))


