import os

class DockerConfig(object):
    # these should be declared in docker-environment.env
    MAIL_SERVER = os.environ.get('MAIL_SERVER') or "smtp.gmail.com"
    MAIL_PORT = int(os.environ.get('MAIL_PORT') or 587)
    MAIL_USE_TLS = os.environ.get('MAIL_USE_TLS') is not None
    MAIL_USERNAME = os.environ.get('MAIL_USERNAME')
    MAIL_PASSWORD = os.environ.get('MAIL_PASSWORD')
    SECRET_KEY = os.environ.get('SECRET_KEY') or "dev"
    ADMINS = ['toxproemail@gmail.com']
