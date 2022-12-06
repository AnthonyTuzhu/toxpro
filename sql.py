from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import smtplib
import string
import random
from ciipro_config import CIIProConfig

def passwordGenerator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

def passwordReset(username, temp_password, new_password, User, db):
    user = User.query.filter_by(username=username).first()
    if not user:
        return 'Incorrect Username'
    if user.check_password(temp_password) == False:
        return 'Incorrect Password'
    user.set_password(new_password)
    db.session.commit()
    return 'Password succesfully changed'

def passwordRetrieval(user_email, User, db):
    user = User.query.filter_by(email=user_email).first()
    if not user:
        return "No email"
    temp_password = passwordGenerator()
    user.set_password(temp_password)
    db.session.commit()
    msg = MIMEMultipart('alternative')
    msg['Subject'] = "Password Reset from CIIPro"
    msg['From'] = "{0}@gmail.com".format(CIIProConfig.CIIPRO_EMAIL)
    password = CIIProConfig.CIIPRO_EMAIL_PW
    username = CIIProConfig.CIIPRO_EMAIL
    msg['To'] = str(user_email)
    text = "Your temporary password for CIIPro is %s Please go to www.ciipro-devel.ccib.rutgers.edu/passreset to reset" % str(temp_password)
    html = """\
    <html>
    <head></head>
    <body>
    <p>
    Your temporary password for CIIPro is %s.  Click <a href='ccib.rutgers.edu/passreset'>here</a> to reset your password.
    </p>
    </body>
    </html>
    """ % str(temp_password)
    
    m_text1 = MIMEText(text, 'plain')
    m_text2 = MIMEText(html, 'html')
    
    msg.attach(m_text1)
    msg.attach(m_text2)
    
    s = smtplib.SMTP('smtp.gmail.com:587')
    s.ehlo()
    s.starttls()
    s.login(username, password)
    s.sendmail("{0}@gmail.com".format(CIIProConfig.CIIPRO_EMAIL), str(user_email), msg.as_string())
    s.quit()

def usernameRetrieval(user_email, User, db):
    user = User.query.filter_by(email=user_email).first()
    if not user:
        return "No email"
    ciipro_username = user.username
    msg = MIMEMultipart('alternative')
    msg['Subject'] = "Username Retrieval from CIIPro"
    msg['From'] = "{0}@gmail.com".format(CIIProConfig.CIIPRO_EMAIL)
    password = CIIProConfig.CIIPRO_EMAIL_PW
    username = CIIProConfig.CIIPRO_EMAIL
    msg['To'] = str(user_email)
    text = "Your username for CIIPro is %s" % str(ciipro_username)
    html = """\
    <html>
    <head></head>
    <body>
    <p>
    Your username for CIIPro is %s
    </p>
    </body>
    </html>
    """ % str(ciipro_username)
    
    m_text1 = MIMEText(text, 'plain')
    m_text2 = MIMEText(html, 'html')
    
    msg.attach(m_text1)
    msg.attach(m_text2)
    
    s = smtplib.SMTP('smtp.gmail.com:587')
    s.ehlo()
    s.starttls()
    s.login(username, password)
    s.sendmail("{0}@gmail.com".format(CIIProConfig.CIIPRO_EMAIL), str(user_email), msg.as_string())
    s.quit()
    
    