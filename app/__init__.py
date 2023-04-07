import os

from flask import Flask
from flask_login import LoginManager
from redis import Redis
import rq


def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SQLALCHEMY_DATABASE_URI='sqlite:///' + os.path.join(app.instance_path, 'toxpro.sqlite'),
        SQLALCHEMY_TRACK_MODIFICATIONS=False,

    )

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_object('app.config.DockerConfig')
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
    migrate.init_app(app, db, render_as_batch=True)
    login_manager = LoginManager()
    login_manager.login_view = 'auth.login'
    login_manager.init_app(app)

    @login_manager.user_loader
    def load_user(user_id):
        return User.query.get(int(user_id))

    # register blueprint:  https://flask.palletsprojects.com/en/2.0.x/tutorial/views/
    from . import auth
    app.register_blueprint(auth.bp)

    from . import cheminf
    app.register_blueprint(cheminf.bp)

    # register main outline: https://flask.palletsprojects.com/en/2.0.x/tutorial/blog/
    from . import toxpro
    app.register_blueprint(toxpro.bp)
    app.add_url_rule('/', endpoint='index')

    #app.redis = Redis.from_url('redis://toxpro-redis-1:6379')
    app.redis = Redis('toxpro-redis-1', 6379)
    app.task_queue = rq.Queue('toxpro-tasks', connection=app.redis)

    from .email import mail
    mail.init_app(app)

    from .errors import not_found, internal_error
    # errors
    # TODO: put these in their own module
    app.register_error_handler(404, not_found)
    app.register_error_handler(500, internal_error)

    import logging
    from logging.handlers import SMTPHandler

    if not app.debug:
        if app.config['MAIL_SERVER']:
            auth = None
            if app.config['MAIL_USERNAME'] or app.config['MAIL_PASSWORD']:
                auth = (app.config['MAIL_USERNAME'], app.config['MAIL_PASSWORD'])
            secure = None
            if app.config['MAIL_USE_TLS']:
                secure = ()
            mail_handler = SMTPHandler(
                mailhost=(app.config['MAIL_SERVER'], app.config['MAIL_PORT']),
                fromaddr='no-reply@' + app.config['MAIL_SERVER'],
                toaddrs=app.config['ADMINS'], subject='ToxPro Error',
                credentials=auth, secure=secure)
            mail_handler.setLevel(logging.ERROR)
            app.logger.addHandler(mail_handler)


    return app


# # needed for flask-migate
# application = create_app()
#
# if __name__ == '__main__':
#     application.run()