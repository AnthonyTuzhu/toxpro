"""

Module that instantiates flask-script manager for running stand-alone scripts on the app

Ref link: https://flask-script.readthedocs.io/en/latest/

"""

from routes import app
from flask_script import Manager
from routes import User

manager = Manager(app)

@manager.command
def get_db_users():
    users = User.query.all()
    print("\n".join(map(str, users)))


if __name__ == "__main__":
    manager.run()
