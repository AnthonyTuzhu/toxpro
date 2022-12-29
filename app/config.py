import os

class DevConfig(object):

    MAIL_SERVER = os.environ.get('MAIL_SERVER') or "smtp.gmail.com"
    MAIL_PORT = int(os.environ.get('MAIL_PORT') or 587)
    MAIL_USE_TLS = os.environ.get('MAIL_USE_TLS') is not None
    MAIL_USERNAME = os.environ.get('MAIL_USERNAME')
    MAIL_PASSWORD = os.environ.get('MAIL_PASSWORD')
    ADMINS = ['toxproemail@gmail.com']
    MAIL_USE_TLS = 1,
    MAIL_USERNAME = 'toxproemail'
    MAIL_PASSWORD = 'umnwwonlugcjbkhl'