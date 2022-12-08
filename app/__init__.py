import os

from flask import Flask
from flask_login import LoginManager
def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        SQLALCHEMY_DATABASE_URI = 'sqlite:///' + os.path.join(app.instance_path, 'toxpro.sqlite'),
        SQLALCHEMY_TRACK_MODIFICATIONS=False
    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    # add the db: https://flask.palletsprojects.com/en/2.0.x/tutorial/database/
    from .db_models import db, User, migrate
    db.init_app(app)
    migrate.init_app(app, db)
    login_manager = LoginManager()
    login_manager.login_view = 'auth.login'
    login_manager.init_app(app)

    @login_manager.user_loader
    def load_user(user_id):
        return User.query.get(int(user_id))

    # register blueprint:  https://flask.palletsprojects.com/en/2.0.x/tutorial/views/
    from . import auth
    app.register_blueprint(auth.bp)

    # register main outline: https://flask.palletsprojects.com/en/2.0.x/tutorial/blog/
    from . import toxpro
    app.register_blueprint(toxpro.bp)
    app.add_url_rule('/', endpoint='index')

    # a simple page that says hello
    @app.route('/hello')
    def hello():
        return 'Hello, World!'

    return app


# # needed for flask-migate
# application = create_app()
#
# if __name__ == '__main__':
#     application.run()