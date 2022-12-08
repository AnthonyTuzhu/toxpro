# this module is outlined
# here: https://flask.palletsprojects.com/en/2.0.x/tutorial/database/
# and help from
# https://www.digitalocean.com/community/tutorials/how-to-add-authentication-to-your-app-with-flask-login
# and https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-v-user-logins

from flask_sqlalchemy import SQLAlchemy
from flask_login import UserMixin
from flask_migrate import Migrate
from werkzeug.security import check_password_hash
import click
from flask import current_app, g
from flask.cli import with_appcontext
from datetime import datetime

db = SQLAlchemy()
migrate = Migrate(db)
class User(db.Model, UserMixin):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(64), index=True, unique=True)
    email = db.Column(db.String(120), index=True, unique=True)
    password_hash = db.Column(db.String(128))
    last_login = db.Column(db.DateTime)
    user_created = db.Column(db.DateTime)
    datasets = db.relationship('Dataset', backref='owner', lazy='dynamic')
    def __repr__(self):
        return '<User {}>'.format(self.username)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)


class Dataset(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))
    dataset_name = db.Column(db.String)

    created = db.Column(db.DateTime, index=True, default=datetime.utcnow)
    __table_args__ = (
        db.UniqueConstraint('user_id', 'dataset_name', name='_user_dataset_uc'),
    )

class Chemical(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    inchi = db.Column(db.String)

# class Chemical(db.Model):
#     p


@click.command('init-db')
@with_appcontext
def init_db_command():
    """Clear the existing data and create new tables."""
    from . import create_app
    db.create_all(app=create_app())
    click.echo('Initialized the database.')