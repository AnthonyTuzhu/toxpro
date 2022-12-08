from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, jsonify, send_from_directory
)
import plotly
import plotly.express as px
import json, os
# from werkzeug.exceptions import abort
#
# from app.auth import login_required
# from app.db import get_db

from app.toxpro_db import DB as toxprodb
from flask_login import current_user, login_required

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
    print("okay")
    return render_template('toxpro/home.html')


@bp.route('/upload', methods=['GET'])
@login_required
def upload():
    """
    displays the homepage

    """

    return render_template('toxpro/upload.html')


@bp.route('/curator', methods=['GET'])
@login_required
def curator():
    """
    displays the homepage

    """

    return render_template('toxpro/curator.html')


@bp.route('/profile', methods=['GET'])
@login_required
def profile():
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