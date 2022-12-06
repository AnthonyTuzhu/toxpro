# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 14:45:16 2015

@author: Daniel
"""
from routes import db

class User(db.Model):
    
    
    __tablename__ = "users"
    id = db.Column('user_id', db.Integer, primary_key=True)
    username = db.Column('username', db.String(20), unique=True, index=True)
    password = db.Column('password', db.String(10))
    email = db.Column('email', db.String(50), unique=True, index=True)
    
    def __init__(self, username, password, email):
        self.username = username
        self.password = password
        self.email = password
        self.email = email
        
    def is_authenticated(self):
        return True
        
    def is_active(self):
        return True
    
    def is_anonymous(self):
        return False
        
    def get_id(self):
        return unicode(self.id)
        
    def __repr__(self):
        return '<User %r>' % (self.username)