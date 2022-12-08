# this module is outlined
# here: https://flask.palletsprojects.com/en/2.0.x/tutorial/views/


import functools
import datetime
from datetime import timezone

from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for
)
from werkzeug.security import check_password_hash, generate_password_hash

from app.db import get_db

bp = Blueprint('auth', __name__, url_prefix='/auth')

@bp.route('/register', methods=('GET', 'POST'))
def register():
    """ Registers a new user.  checkRecaptcha() must return True to register user.


    """
    if request.method == 'GET':
        return render_template('auth/register.html')
    print("test")
    username = request.form['username'].strip()
    password = request.form['password'].strip()
    email = request.form['email'].strip()
    error = None

    db = get_db()
    if not username:
        error = 'Username is required.'
    elif not password:
        error = 'Password is required.'

    # check to make sure password if valid
    password_check = PasswordCheck(password_string=password)
    if not password_check.has_numbers():
        error = 'Password must contain at least one number.'
    elif not password_check.has_letters():
        error = 'Password must contain at least one letter.'
    elif not password_check.is_n_letters_long(n=5):
        error = 'Password must be at least 5 characters long.'

    # check to make sure email is valid
    # TODO: add step to confirm registration
    email_check = EmailCheck(email_string=email)
    if not email_check.is_valid():
        error = 'Email is not valid.'


    if error is None:
        try:
            dt = datetime.datetime.now(timezone.utc)

            insert_values = (username,
                             generate_password_hash(password),
                             email,
                             dt,
                             dt)
            db.execute(
                "INSERT INTO user (username, password, email, created_date, last_login) VALUES (?, ?, ?, ?, ?)",
                insert_values,
            )
            db.commit()
        except db.IntegrityError:
            user = db.execute(
                'SELECT * FROM user WHERE username = ?', (username,)
            ).fetchone()

            email_result = db.execute(
                'SELECT * FROM user WHERE email = ?', (email,)
            ).fetchone()

            if user:
                error = f"User {username} is already registered."
            elif email_result:
                error = f"Email ({email}) is already registered."

            if (user and email_result):
                error = 'Email and username already registered.'

        else:
            flash('Successfully registered user', 'success')
            return redirect(url_for("auth.login"))

    flash(error, 'danger')
    return redirect(url_for('auth.register'))


@bp.route('/login', methods=('GET', 'POST'))
def login():
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        db = get_db()
        error = None
        user = db.execute(
            'SELECT * FROM user WHERE username = ?', (username,)
        ).fetchone()

        if user is None:
            error = 'Incorrect username.'
        elif not check_password_hash(user['password'], password):
            error = 'Incorrect password.'


        if error is None:
            session.clear()
            session['user_id'] = user['id']
            dt = datetime.datetime.now(timezone.utc)
            sql = ''' UPDATE user 
                      SET last_login = ? 
                      WHERE id = ? '''
            db.execute(sql, (dt, user['id']))
            db.commit()
            return redirect(url_for('index'))

        flash(error)

    return render_template('auth/login.html')

@bp.before_app_request
def load_logged_in_user():
    user_id = session.get('user_id')

    if user_id is None:
        g.user = None
    else:
        g.user = get_db().execute(
            'SELECT * FROM user WHERE id = ?', (user_id,)
        ).fetchone()

@bp.route('/logout')
def logout():
    session.clear()
    return redirect(url_for('index'))


def login_required(view):
    @functools.wraps(view)
    def wrapped_view(**kwargs):
        if g.user is None:
            return redirect(url_for('auth.login'))

        return view(**kwargs)

    return wrapped_view


class PasswordCheck:
    """ class to check password criteria """

    def __init__(self, password_string: str):
        self.password_string = password_string

    def is_n_letters_long(self, n: int):
        return len(self.password_string) >= n

    def has_numbers(self):
        return any(letter.isdigit() for letter in self.password_string)

    def has_letters(self):
        return not all(letter.isdigit() for letter in self.password_string)

class EmailCheck:

    """ class to check email criteria """

    def __init__(self, email_string: str):
        self.email_string = email_string

    def is_valid(self):
        has_at = '@' in self.email_string
        return has_at