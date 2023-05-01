# This is the module for mapping database ORM models
# to the SQLite database
#
# this module is outlined
# here: https://flask.palletsprojects.com/en/2.0.x/tutorial/database/
# and help from
# https://www.digitalocean.com/community/tutorials/how-to-add-authentication-to-your-app-with-flask-login
# and https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-v-user-logins
#
# There was also a nice tutorial on how register a user
# outlined here: https://realpython.com/handling-email-confirmation-in-flask/

from flask_sqlalchemy import SQLAlchemy
from flask_login import UserMixin
from flask_migrate import Migrate
from werkzeug.security import check_password_hash, generate_password_hash
import click
from flask import current_app, g
from flask.cli import with_appcontext
from datetime import datetime
from sqlalchemy import MetaData
import redis
import rq
import jwt
from time import time


from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

# from app import create_app

# app = create_app()
# app.app_context().push()

# this is a fix for unique constrains
# and forien keys in Flask-migrate
# see SO: https://stackoverflow.com/questions/45527323/flask-sqlalchemy-upgrade-failing-after-updating-models-need-an-explanation-on-h

naming_convention = {
    "ix": 'ix_%(column_0_label)s',
    "uq": "uq_%(table_name)s_%(column_0_name)s",
    "ck": "ck_%(table_name)s_%(column_0_name)s",
    "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
    "pk": "pk_%(table_name)s"
}
db = SQLAlchemy(metadata=MetaData(naming_convention=naming_convention))

migrate = Migrate(db)

class User(db.Model, UserMixin):
    """ Main class to handle users """
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(64), index=True, unique=True)
    email = db.Column(db.String(120), index=True, unique=True)
    password_hash = db.Column(db.String(128))
    last_login = db.Column(db.DateTime)
    user_created = db.Column(db.DateTime)

    # an attribute to restrict functionality
    admin = db.Column(db.Boolean, default=False)

    confirmed = db.Column(db.Boolean, default=False)
    confirmed_on = db.Column(db.DateTime)

    datasets = db.relationship('Dataset', backref='owner', lazy='dynamic', cascade="all, delete-orphan")
    tasks = db.relationship('Task', backref='user', lazy='dynamic', cascade="all, delete-orphan")
    def __repr__(self):
        return '<User {}>'.format(self.username)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

    def launch_task(self, name, description, *args, **kwargs):
        rq_job = current_app.task_queue.enqueue('app.tasks.' + name,
                                                *args, **kwargs)
        task = Task(id=rq_job.get_id(), name=name, description=description,
                    user=self)
        db.session.add(task)
        return task

    def get_tasks_in_progress(self):
        return Task.query.filter_by(user=self, complete=False).all()

    def get_task_in_progress(self, name):
        return Task.query.filter_by(name=name, user=self,
                                    complete=False).first()

    def get_token(self, kind: str='reset_password', expires_in=600):
        return jwt.encode(
            {kind: self.id, 'exp': time() + expires_in},
            current_app.config['SECRET_KEY'], algorithm='HS256')

    @staticmethod
    def verify_token(token, kind):
        try:
            id = jwt.decode(token, current_app.config['SECRET_KEY'],
                            algorithms=['HS256'])[kind]
        except:
            return
        return User.query.get(id)

class Dataset(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id', name='user-dataset'))
    dataset_name = db.Column(db.String)
    chemicals = db.relationship('Chemical', backref='dataset', lazy='dynamic', cascade="all, delete-orphan")
    qsar_models = db.relationship('QSARModel', backref='dataset', lazy='dynamic', cascade="all, delete-orphan")

    created = db.Column(db.DateTime, default=datetime.utcnow)
    __table_args__ = (
        db.UniqueConstraint('user_id', 'dataset_name', name='_user_dataset_uc'),
    )

    def get_chemicals(self):
        return Chemical.query.join(Dataset,
                              Dataset.id == Chemical.dataset_id).filter(Dataset.id == self.id)

class Chemical(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    inchi = db.Column(db.String)
    dataset_id = db.Column(db.Integer, db.ForeignKey('dataset.id'))
    activity = db.Column(db.Integer)
    compound_id = db.Column(db.String)

    def to_dict(self, structure_as_svg=False):
        result = {
            'Chemical': self.compound_id,
            'Activity': self.activity,
            'Structure': self.inchi
        }

        if structure_as_svg:
            result['Structure'] = self.get_svg()

        return result

    def get_svg(self):
        mol = Chem.MolFromInchi(self.inchi)
        d2d = rdMolDraw2D.MolDraw2DSVG(current_app.config['SVG_DISPLAY_WIDTH'],
                                       current_app.config['SVG_DISPLAY_HEIGHT'])
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        return d2d.GetDrawingText()


class QSARModel(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id', name='user-dataset'))
    name = db.Column(db.String)
    algorithm = db.Column(db.String)
    descriptors = db.Column(db.String)
    dataset_id = db.Column(db.Integer, db.ForeignKey('dataset.id'))
    sklearn_model = db.Column(db.BLOB)
    created = db.Column(db.DateTime, default=datetime.utcnow)

    cvresults = db.relationship('CVResults',
                                backref='qsar_model',
                                cascade="all, delete-orphan",
                                uselist=False # specifies one-to-one
                                )

class CVResults(db.Model):
    """ class for handling and storing five fold cross validation results from QSAR modeling

    child of QSARModel

     """
    id = db.Column(db.Integer, primary_key=True)
    accuracy = db.Column(db.Float)
    f1_score = db.Column(db.Float)
    area_under_roc = db.Column(db.Float)
    cohens_kappa = db.Column(db.Float)
    michaels_correlation = db.Column(db.Float)
    precision = db.Column(db.Float)
    recall = db.Column(db.Float)
    specificity = db.Column(db.Float)
    correct_classification_rate = db.Column(db.Float)

    qsar_model_id = db.Column(db.Integer, db.ForeignKey('qsar_model.id'))


class Task(db.Model):
    id = db.Column(db.String(36), primary_key=True)
    name = db.Column(db.String(128), index=True)
    description = db.Column(db.String(128))
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))
    complete = db.Column(db.Boolean, default=False)
    time_submitted = db.Column(db.DateTime, default=datetime.utcnow)
    time_completed = db.Column(db.DateTime)

    def get_rq_job(self):
        try:
            rq_job = rq.job.Job.fetch(self.id, connection=current_app.redis)
        except (redis.exceptions.RedisError, rq.exceptions.NoSuchJobError):
            return None
        return rq_job

    def get_progress(self):
        job = self.get_rq_job()
        return job.meta.get('progress', 'Queued') if job is not None else 'Finished'

def create_db(overwrite=False):
    """Clear the existing data and create new tables."""
    from . import create_app
    import os
    app = create_app()
    if overwrite and os.path.exists(os.path.join(app.instance_path, 'toxpro.sqlite')):
        os.remove(os.path.join(app.instance_path, 'toxpro.sqlite'))
    db.create_all(app=app)

@click.command('init-db')
@with_appcontext
def init_db_command():
    """Clear the existing data and create new tables."""
    from . import create_app
    db.create_all(app=create_app())
    click.echo('Initialized the database.')