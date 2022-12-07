# -*- coding: utf-8 -*-
import os, sys

filename = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

sys.path.insert(0, filename)
print(filename)


import os
import ntpath

from flask import Flask, render_template, flash, session, \
    redirect, url_for, g, send_file, request, jsonify
from flask_login import login_user, logout_user, current_user, login_required, LoginManager
from werkzeug.utils import secure_filename
from werkzeug.security import generate_password_hash, check_password_hash
from flask_sqlalchemy import SQLAlchemy
from sql import passwordRetrieval, usernameRetrieval, passwordReset

import json, glob

#import urllib
import zipfile

import sqlalchemy.ext


import pandas as pd
import numpy as np


from ml import get_class_stats



# These variables are configured in CIIProConfig
# TODO: put all this in a true config file
app = Flask(__name__)
#app.config['UPLOAD_FOLDER'] = CIIProConfig.UPLOAD_FOLDER
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.db'
app.secret_key = 'test'
app.config['RECAPTCHA_PRIVATE_KEY'] = 'test'

# register the routes api

db = SQLAlchemy(app)
login_manager = LoginManager()
login_manager.init_app(app)
login_manager.login_view = 'login'


class User(db.Model):
    tablename__ = "users"
    id = db.Column('user_id', db.Integer, primary_key=True)
    username = db.Column('username', db.String(20), unique=True, index=True)
    pw_hash = db.Column('password', db.String(10))
    email = db.Column('email', db.String(50), unique=True, index=True)
    
    def __init__(self, username, password, email):
        self.username = username
        self.set_password(password)
        self.email = email
        
    def set_password(self, password):
        self.pw_hash = generate_password_hash(password, method='pbkdf2:sha1', 
                                              salt_length=8)
        
    def check_password(self, password):
        return check_password_hash(self.pw_hash, password)
        
    def is_authenticated(self):
        return True
        
    def is_active(self):
        return True
    
    def is_anonymous(self):
        return False
        
    def get_id(self):
        return str(self.id)
        
    def __repr__(self):
        return '<User %r>' % (self.username)

    # I added these functions,
    # very helpful to probably bring a lot of stuff into this class

    def get_user_folder(self, folder_name):

        # make folder if it doesnt exists
        folder = os.path.join(CIIProConfig.UPLOAD_FOLDER, self.username, folder_name)
        if not os.path.exists(folder):
            os.makedirs(folder)
        return folder

    def get_user_datasets(self, set_type="training"):
        """ returns the datasets for a users """
        return ds_io.get_datasets_names_for_user(self.get_user_folder('datasets'), set_type=set_type)

    def get_user_bioprofiles(self):
        """ returns profiles for a user """
        return ciipro_io.get_profiles_names_for_user(self.get_user_folder('profiles'))

    def get_user_fp_profiles(self):
        """ returns profiles for a user """
        return ciipro_io.get_profiles_names_for_user(self.get_user_folder('fp_profiles'))

    def load_dataset(self, ds_name):

        """ load a dataset object given a dataset name for a user """

        if (ds_name in self.get_user_datasets("training")) or (ds_name in self.get_user_datasets("test")):
            ds_json_file = os.path.join(self.get_user_folder('datasets'), '{}.json'.format(ds_name))

            return ds.DataSet.from_json(ds_json_file)


    def load_bioprofile(self, bp_name):

        """ load a dataset object given a dataset name for a user """

        if bp_name in self.get_user_bioprofiles():
            bp_json_file = os.path.join(self.get_user_folder('profiles'), '{}.json'.format(bp_name))

            return bp.Bioprofile.from_json(bp_json_file)

    def load_adj_matrix(self, clustering_name):
        adj_json_file = os.path.join(self.get_user_folder('fp_profiles'),
                                     '{}_adj_matrix.json'.format(clustering_name))
        return fp.AdjMatrix.from_json(adj_json_file)

    def load_fp_profile(self, clustering_name):
        fp_json_file = os.path.join(self.get_user_folder('fp_profiles'),
                                     '{}.json'.format(clustering_name))
        return fp.FPprofile.from_json(fp_json_file)


    
#db.create_all()
@login_manager.user_loader
def load_user(id):
    return User.query.get(int(id))

@app.after_request
def add_header(response):
    """
    Add headers to both force latest IE rendering engine or Chrome Frame,
    and also to cache the rendered page for 10 minutes.
    """
    response.headers['X-UA-Compatible'] = 'IE=Edge,chrome=1'
    response.headers['Cache-Control'] = 'public, max-age=0'
    return response


@app.route('/') 
def home():
    return render_template('home.html')

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'GET':
        return render_template('login.html')
        
    username = request.form['username'].strip()
    password = request.form['password'].strip()
    registered_user = User.query.filter_by(username=username).first()
    
    if registered_user is None:
        flash('Username does not exist', 'danger')
        return redirect(url_for('login'))
       
    if registered_user.check_password(password) == False:
        flash('Password does not match user', 'danger')
        return redirect(url_for(''))
    
    login_user(registered_user)
    flash('Logged in successfully', 'info')
    return redirect(request.args.get('next') or url_for('home'))

@app.route('/logout')
@login_required
def logout():
    """ Logs out user and returns them to the homepage.
    
    """
    
    session.pop('compound_file', None)
    logout_user()
    flash('Logged out successfully', 'info')
    return redirect(url_for('home'))

def checkRecaptcha(response, secretkey):
    """ Checks the response to the recaptcha entry on the login page. Returns True if response == recaptcha diplayed
        on the website.
    
    response (str): User supplied recaptcha key.
    secretkey (str): Secret key for the site supplied by Google.
    """
    url = 'https://www.google.com/recaptcha/api/siteverify?'      
    url = url + 'secret=' + secretkey
    url = url + '&response=' + response
    try:
        jsonobj = json.loads(urllib.request.urlopen(url).read().decode('utf-8'))
        if jsonobj['success']:
            return True
        else:
            return False
    except Exception as e:
        print(e)

        return False
      
      
@app.route('/register', methods=['GET', 'POST']) 
def register():
    """ Registers a new user.  checkRecaptcha() must return True to register user.
        Upon successful registration, creates root directory. And four subdirectory folders:
            
            Username
               ├───compounds
               ├───biosims
               ├───profiles
               └───converter   
    """
    if request.method == 'GET':
        return render_template('register.html')

    try:
        user = User(request.form['username'], request.form['password'], request.form['email'])
    except sqlalchemy.exc.IntegrityError as e:
        flash('Sorry, either user exists with that username or email.', 'danger')
    recaptcha = request.form['g-recaptcha-response']
    if checkRecaptcha(recaptcha, CIIProConfig.SECRET_KEY):
        db.session.add(user)
        db.session.commit()
        directory = CIIProConfig.UPLOAD_FOLDER + '/' + str(user.username)
        if not os.path.exists(directory):
            os.makedirs(directory)
            sub_folders = ["compounds", "biosims", "profiles", "converter", "test_sets", "NNs", "datasets"]
            for sf in sub_folders:
                sub_directory = directory + '/' + sf
                os.makedirs(sub_directory)
        flash('User successfully registered')
    else:
        flash('Registration failed', 'danger')
    return redirect(url_for('login'))

@app.route('/passwordrecovery', methods=['GET', 'POST']) 
def passrecov():
    """ Checks to see if an email is associated with a user to recover password.  If so, returns emailsent.html.
  
    """
    if request.method == 'GET':
        return render_template('passwordrecovery.html')
    
    if request.method == 'POST':
        email = request.form['email']
        response = passwordRetrieval(email, User, db)
        if response == "No email":
            error = "No user associated with that email, please register"
            return render_template('passwordrecovery.html')
        else:
            return render_template('emailsent.html')

@app.route('/usernamerecovery', methods=['GET', 'POST']) 
def usernamerecov():
    """ Checks to see if an email is associated with a user to recover username.  If so, returns emailsent.html.
  
    """
    if request.method == 'GET':
        return render_template('usernamerecovery.html')
    if request.method == 'POST':
        email = request.form['email']
        response = usernameRetrieval(email, User, db)
        if response == "No email":
            error = "No user associated with that email, please register"
            return render_template('passwordrecovery.html', error=error)
        else:
            return render_template('emailsent.html')            

@app.route('/passreset', methods=['GET', 'POST'])
def passreset():
    """ Resets a users password.
    
    """
    if request.method == 'GET':
        return render_template('passreset.html')
    
    if request.method == 'POST':
        username = request.form['username']
        temp_password = request.form['temp_password']
        new_password = request.form['new_password']
        conf_password = request.form['conf_password']
        
        if new_password != conf_password:
            flash('New Password and confirmation don\'t match', 'danger')
            return render_template('passreset.html')
        
        response = passwordReset(username, temp_password, new_password, User, db)
        if response == 'Password succesfully changed':
            registered_user = User.query.filter_by(username=username).first()
            login_user(registered_user)
            flash(response, 'info')
            return render_template('home.html')
        
        return render_template('passreset.html')
   
        
@app.before_request
def before_request():
    g.user = current_user

@app.route('/tutorial')
def tutorial():
    return render_template('tutorial.html')

@app.route('/sendbiopro')
def sendbiopro():
      return render_template('StatsGlossary.html')

@app.route('/sendbiosim')
def sendbiosim():
      return send_file('tutorial_samples/BioSim.txt', as_attachment=True)      

@app.route('/sendbioconf')
def sendbiosimconf():
      return send_file('tutorial_samples/BioSim_Conf.txt', as_attachment=True)      

@app.route('/sendbionn')
def sendbionn():
      return send_file('tutorial_samples/Bioneighbor.txt', as_attachment=True)
      
@app.route('/sendactivity')
def sendactivity():
      return send_file('tutorial_samples/Activity.txt', as_attachment=True)
      
@app.route('/sendbiopred')
def sendbiopred():
      return send_file('tutorial_samples/BioPred.txt', as_attachment=True)


@app.route('/statsglossary')
def statsglossary():
    """ Displays statsglossary page. 
    
    """    
    return render_template('statsglossary.html')
    
@app.route('/datasets', methods=['GET', 'POST'])
@login_required
def datasets():
    """ Displays datasets page with all available datasets in users compound folder. 
    
    """
    datasets = g.user.get_user_datasets(set_type='training') + g.user.get_user_datasets(set_type='test')

    return render_template('datasets.html', datasets=datasets,
                           testsets=g.user.get_user_datasets(set_type='test'),
                          username=g.user.username)
                           
                           
@app.route('/uploaddataset', methods=['POST', 'GET'])
@login_required
def uploaddataset():
    """ Uploads a file from user and saves to users' compounds folder.  Converts non-PubChem CID identifiers to CIDS using 
        PubChem's PUG Rest.
        
        Requests:
            input_type (str): radio button from page
            compound_file: user file upload, first column should be compounds, second should be activity.
            model_type: training or test set upload
    """
   

    
    username = g.user.username

    # requests
    input_type = request.form['input_type'].lower()
    file = request.files['compound_file']
    model_type = request.form['model_type']
    
    if file and allowed_file(file.filename):
        compound_filename = secure_filename(file.filename)

        user_datasets_folder = g.user.get_user_folder('datasets')
        user_uploaded_file = os.path.join(user_datasets_folder, compound_filename)

        # make the name of the dataset just the basename of the file

        name = ntpath.basename(user_uploaded_file).split('.')[0]

        # I think we have to save this in order to use it, not sure if we car read it otherwise
        file.save(user_uploaded_file)
        print(user_uploaded_file)
        # TODO: Need to write some checks to make sure everything in the uploaded dataset is good
        # TODO: Things like all activies are there, columns match, etc
        # TODO: the following type conversions are necessary for JSON serialization

        identifiers, activities = ciipro_io.parse_upload_file(user_uploaded_file)

        activities = list(map(int, activities))

        if input_type == 'cid':
            identifiers = list(map(int, identifiers))
        else:
            identifiers = list(map(str, identifiers))

        ds_io.write_ds_to_json(identifiers, activities, user_datasets_folder, name, input_type, set_type=model_type)
        os.remove(user_uploaded_file)

        return redirect(url_for('datasets'))  
    else:
        error = "Please attach file"
        return render_template('datasets.html',
                               datasets=g.user.get_user_datasets(set_type='training'),
                               testsets=g.user.get_user_datasets(set_type='test'),
                               username=username, error=error)
     

@app.route('/deletetestset', methods=['POST'])
@login_required
def deletetestset():
    """ Deletes a test set from a users' test sets folder.  
    
        Requests:
            testset_filename (str): radiobutton from datasets page.  
    """
    testset_filename = '{}.json'.format(str(request.form['testset_filename']))
    os.remove(os.path.join(g.user.get_user_folder('datasets'), testset_filename))
    return redirect(url_for('datasets'))  
                            
@app.route('/deletedataset', methods=['POST'])
@login_required
def deletedataset():
    """ Deletes a dataset from a users' compounds folder.  
    
        Requests:
            compound_filename (str): radiobutton from datasets page.  
    """
    compound_filename = '{}.json'.format(str(request.form['compound_filename']))
    os.remove(os.path.join(g.user.get_user_folder('datasets'), compound_filename))
    return redirect(url_for('datasets'))

@app.route('/CIIProfiler') 
@login_required
def CIIProfiler():
    """ Displays CIIProfiler page with all available datasets in users compound folder.
    
    """
    return render_template('CIIProfiler.html', username=g.user.username,
                           datasets=g.user.get_user_datasets(set_type='training'))


@app.route('/CIIProdatasets')
def CIIProdatasets():
    """ Displays CIIProfiler page with all available datasets in users compound folder.

    """
    return render_template('CIIProdatasets.html')


@app.route('/CIIPPredictor')
@login_required
def CIIPPredictor():
    """ Displays CIIPBioNN page with all available profiles in users' profile folder.
    
    """

    return render_template('CIIPPredictor.html', profiles=g.user.get_user_bioprofiles(),
                           username=g.user.username, testsets=g.user.get_user_datasets(set_type='test'))	


# TODO Create a save button to save optimized profile, allow for name
@app.route('/CIIPro_Optimizer', methods=['POST', 'GET'])
@login_required
def CIIPro_Optimizer():
    """ Renderts optimize bioprofile
    """


    return render_template('CIIPro_Optimizer.html',
                           profiles=g.user.get_user_bioprofiles(),
                           username=g.user.username)

def allowed_file(filename): #method that checks to see if upload file is allowed
    return '.' in filename and filename.rsplit('.', 1)[1] in CIIProConfig.ALLOWED_EXTENSIONS


@app.route('/CIIPro_Cluster', methods=['GET', 'POST'])
@login_required
def CIIPro_Cluster():
    """
    """

    if request.method == 'GET':
        return render_template('CIIProCluster.html', profiles=g.user.get_user_bioprofiles(),
                                                    clusters=g.user.get_user_fp_profiles())
    else:

        profile_filename = request.form['profile_filename'].strip()
        clustering_filename = request.form['clustering_filename'].strip()

        threshold = 0.05

        training_profile = g.user.load_bioprofile(profile_filename)

        # load the associated training set
        training_set = g.user.load_dataset(training_profile.meta['training_set'])

        bioprofile = training_profile.to_frame()
        fps = training_set.get_pubchem_fps().loc[bioprofile.index]

        correlations = in_vitro_fingerprint_correlations(bioprofile, fps, threshold=threshold, binarize=True)

        fp_profile = fp.FPprofile.from_dict(correlations)
        fp_profile.name = clustering_filename


        meta = {
            'profile_used': training_profile.name,
            'threshold': threshold
        }

        fp_profile.meta = meta

        fp_profile.to_json(g.user.get_user_folder('fp_profiles'))

        adj_matrix = fp_profile.get_adjacency()

        adj_matrix.to_json(g.user.get_user_folder('fp_profiles'))

        return render_template('CIIProCluster.html', profiles=g.user.get_user_bioprofiles(),
                               clusters=g.user.get_user_fp_profiles())

@app.route('/optimizeassays', methods=['POST'])
@login_required
def removeassays():
    remove_or_keep = request.form['Remove_or_Keep']
    assays = list(map(int, request.form.getlist('aids')))
    profile = bioprofile_to_pandas(session['cur_prof_dir'])
    if remove_or_keep == 'Remove':
        profile.drop(assays, axis=1, inplace=True)
    else:
        profile = profile[assays]

    new_name = session['cur_prof_dir'].replace('.txt',
                                               '_optimized.txt')

    profile.to_csv(new_name, sep='\t')
    profile.to_csv(new_name.replace('biosims', 'profiles'), sep='\t')

    session['cur_prof_dir'] = new_name
    return redirect(url_for('optimizeprofile'))

@app.route('/ciiprofile',  methods=['POST'])
@login_required
def CIIProfile():
    """ Creates a bioprofile from a users selected dataset.  If the bioprofile is not too large, will create a 
        heatmap to display on the website.  
    
        Requests:
            compound_filename (str): Selected dataset from radio button.  
            noOfActives (str): String number of minumum actives per assays from pull down.
            profile_filename (str): Name of profile that will be generated. From text field. 
            sort_by: Column to sort in vitro, in vivo correlations by
                
    """

    if request.method == 'POST':

        ds_name = str(request.form['compound_filename'])
        ds = g.user.load_dataset(ds_name)
        min_actives = int(request.form['noOfActives'])
        profile_name = str(request.form['profile_filename']).strip()

        bioprofile_json = ds.get_bioprofile(min_actives=min_actives)

        act = ds.get_activities(use_cids=True)
        # TODO a better way to tdo this

        bioprofile = bp.Bioprofile(profile_name,
                                   bioprofile_json['cids'],
                                   bioprofile_json['aids'],
                                   bioprofile_json['outcomes'], None, None)

        profile_matrix = bioprofile.to_frame()

        stats_df = getIVIC(act, profile_matrix)
        stats_df.reset_index(inplace=True)
        stats_df = stats_df.rename(str, columns={"index": "aid"})

        meta = {}
        meta['training_set'] = str(ds_name)
        meta['num_total_actives'] = int((profile_matrix == 1).sum().sum())
        meta['num_total_inactives'] = int((profile_matrix == -1).sum().sum())
        meta['num_cmps'] = int(profile_matrix.shape[0])
        meta['num_aids'] = int(profile_matrix.shape[1])


        # this has to happen to store AIDs as attributes in JSON
        stats_df.index = list(map(str, stats_df.index))

        bioprofile = bp.Bioprofile(profile_name,
                                   bioprofile_json['cids'],
                                   bioprofile_json['aids'],
                                   bioprofile_json['outcomes'], stats_df.to_dict('records'), meta)


        bioprofile.to_json(g.user.get_user_folder('profiles'))



        profile_matrix = bioprofile.to_frame()


        flash('Success! A profile was created consisting '
              'of {0} compounds and {1} bioassays'.format(profile_matrix.shape[0], profile_matrix.shape[1]), 'info')
        return render_template('CIIProfiler.html', stats=None, datasets=g.user.get_user_datasets(set_type='training'))


    
@app.route('/CIIPPredict', methods=['POST'])
@login_required
def CIIPPredict():
    """ Creates a biosimilarity and biological nearest neighbors from a user selected profile. 
        If an activity file is uploaded, will calculate and display in vitro, in vivo correlations as well as the 
        results of leave one out predictions of the test set.  
        
        Form Requests: 
            profile_filename: The name of the profile to use for biosimilarity and biological NNs.  From radio buttons.
            cutoff: Biosimilarity cutoff to use for NN calculation.  From text field.
            confidence: Minimum confidence score to use for NN calculation.  From text field.
            submit: Type of submit button, either Delete or Submit
    """

    if request.method == 'POST':


        profile_filename = request.form['profile_filename']
        compound_filename = request.form['compound_filename']

        biosim_cutoff = float(request.form['cutoff'])
        conf_cutoff = float(request.form['conf'])
        nn_cutoff = float(request.form['nns'])

        training_profile = g.user.load_bioprofile(profile_filename)

        training_matrix = training_profile.to_frame()

        training_data = g.user.load_dataset(training_profile.meta['training_set'])

        test_data = g.user.load_dataset(compound_filename)



        # get a full test set profile

        test_profile_json = test_data.get_bioprofile()

        # I guess since a lot of these dont get used might be better to
        # create two classes one for training profiles and one for test
        test_profile = bp.Bioprofile(None,
                                     test_profile_json['cids'],
                                     test_profile_json['aids'],
                                     test_profile_json['outcomes'],
                                     None,
                                     None)

        test_matrix = test_profile.to_frame()

        # only use the intersection of the test profile and the training profile
        shared_assays = training_matrix.columns.intersection(test_matrix.columns)

        # TODO: Need to come up with an error if no test compounds have assays in the training profile

        # align both axis
        test_matrix = test_matrix.loc[:, shared_assays]
        training_matrix = training_matrix.loc[:, shared_assays]



        trainin_activites = training_data.get_activities(use_cids=True)
        trainin_activites = trainin_activites[training_matrix.index]

        # currently weight is ratio of actives:inactives in training profile
        # but never have it above 1
        # TODO: add this as a paramater to let users choose

        act_to_inact_ratio =(training_matrix == 1).sum().sum() / (training_matrix == -1).sum().sum()

        inact_weight = act_to_inact_ratio if act_to_inact_ratio <= 1 else 1

        biodis, conf = biosim.biosimilarity_distances(test_matrix.values,
                                                      training_matrix.values,
                                                      weight=inact_weight)

        biosimilarity = 1-biodis

        bio_nns = biosim.get_k_bioneighbors(biosimilarity, conf,
                                        k=nn_cutoff,
                                        biosim_cutoff=biosim_cutoff,
                                        conf_cutoff=conf_cutoff)

        train_fps = training_data.get_pubchem_fps()
        test_fps = test_data.get_pubchem_fps()

        # make sure axis are aligned
        train_fps = train_fps.loc[training_matrix.index]
        test_fps = test_fps.loc[test_matrix.index]


        # now get chem similarity
        chemdis = biosim.chem_distances(test_fps.values, train_fps.values)
        chemsim = 1-chemdis
        chem_nns = biosim.get_k_chem_neighbors(chemsim, k=nn_cutoff)

        # now package is all up in a json ob

        data = []

        for i, test_cmp in enumerate(test_matrix.index):

            test_cmp_data = {}
            test_cmp_data["cid"] = int(test_cmp)
            test_cmp_data['bionn'] = list(map(int, training_matrix.index[bio_nns[i]]))
            test_cmp_data['bioSims'] = list(map(float, biosimilarity[i, bio_nns[i]]))
            test_cmp_data['bioConf'] = list(map(float, conf[i, bio_nns[i]]))
            test_cmp_data['bioPred'] = float(trainin_activites[training_matrix.index[bio_nns[i]]].mean())

            test_cmp_data['chemnn'] = list(map(int, training_matrix.index[chem_nns[i]]))
            test_cmp_data['chemSims'] = list(map(float, chemsim[i, chem_nns[i]]))
            test_cmp_data['chemPred'] = float(trainin_activites[training_matrix.index[chem_nns[i]]].mean())
            data.append(test_cmp_data)


        biopreds = pd.Series([data_dic["bioPred"] for data_dic in data],
                             index=[data_dic["cid"] for data_dic in data])

        chempreds = pd.Series([data_dic["chemPred"] for data_dic in data],
                             index=[data_dic["cid"] for data_dic in data])


        hybridpreds = (biopreds + chempreds) / 2
        hybridpreds = hybridpreds.apply(np.ceil)

        y_test = test_data.get_activities(use_cids=True).loc[biopreds.index]

        chem_stats = [{"metric": metric, "value": round(value, 2)}
                      for metric, value in get_class_stats(chempreds, y_test).items()]
        bio_stats = [{"metric": metric, "value": round(value, 2)}
                      for metric, value in get_class_stats(biopreds, y_test).items()]
        hybrid_stats = [{"metric": metric, "value": round(value, 2)}
                      for metric, value in get_class_stats(hybridpreds, y_test).items()]


        stats = [bio_stats, chem_stats, hybrid_stats]
        print(stats)
        return render_template('CIIPPredictor.html', profiles=g.user.get_user_bioprofiles(),
                               testsets=g.user.get_user_datasets(set_type='test'),
                               data=data, stats=stats)



@app.route('/similarity<cid>', methods=['GET', 'POST'])
@login_required
def similarity(cid):
    USER_FOLDER = CIIProConfig.UPLOAD_FOLDER + '/' + g.user.username
    USER_TEST_SET_FOLDER = USER_FOLDER +'/test_sets'
    USER_NN_FOLDER = USER_FOLDER + '/NNs'
    
    df = nn_to_pandas(USER_NN_FOLDER + '/' + session['test_set'] + '/' + str(cid) + '.csv')
    #sim_graph = createSimilarityGraph(int(cid), df,  int(session['nns']))
    sim_graph_pic = sim_graph(int(cid), df, int(session['nns']), int(session['max_conf']))
    #cids = session['cids']
    #preds = session['preds']
    #acts = session['acts']
    #len_cids = len(cids)
       
    return render_template('similarity.html', sim_graph=sim_graph_pic)


@app.route('/CIIProTools', methods=['GET', 'POST']) 
@login_required
def CIIProTools():
    """ Displays CIIProTools page with all available training sets
        in users' training folder.
    
    """

    USER_COMPOUND_FOLDER = CIIProConfig.UPLOAD_FOLDER + '/' + g.user.username + '/compounds'
    datasets = [dataset for dataset in os.listdir(USER_COMPOUND_FOLDER)]
    return render_template('CIIProTools.html', datasets=g.user.get_user_datasets(set_type='training'), 
                           username=g.user.username)	

@app.route('/activitycliffs', methods=['GET', 'POST']) 
@login_required
def activitycliffs():
    """ Method Identifies Activity Cliffs in training set
    
    """

    USER_COMPOUND_FOLDER = CIIProConfig.UPLOAD_FOLDER + '/' + g.user.username + '/compounds'
    datasets = [dataset for dataset in os.listdir(USER_COMPOUND_FOLDER)]
    compound_filename = request.form['compound_filename']
    compound_filename = str(compound_filename)
    compound_directory = USER_COMPOUND_FOLDER + '/' + compound_filename
    
    df = activity_cliffs(compound_directory)
    #df.to_csv(compound_directory.replace('compounds', NNs), sep='\t')
    #cliff = cliffTable_bokeh(df)
    df.index.name = 'Target_CID'
        
    writer = pd.ExcelWriter(compound_directory.replace('compounds', 'NNs').replace('.txt', '.xlsx'))
    df.to_excel(writer, sheet_name='Activity Cliffs')
    writer.save()
    session['cur_ciff_dir'] = compound_directory.replace('compounds', 'NNs').replace('.txt', '.xlsx')
        
    return render_template('CIIProTools.html', datasets=g.user.get_user_datasets(set_type='training'), 
                           username=g.user.username, ac=df.to_html())	


@app.route('/sendtutorial')
def trainingsettutorial():
    return send_file('resources/ER_tutorial.zip')
  
 
@app.route('/contact')
def contact():
      return render_template('contact.html')



@app.errorhandler(500)
def internalServiceError(e):
    return render_template('500.html'), 500


# this is where the RESTFul API will be
# not sure if the login_required decorator
# actually works for these

@login_required
@app.route('/get_dataset_overview/<dataset_name>')
def get_dataset_overview(dataset_name):
    print(dataset_name)
    ds = g.user.load_dataset(dataset_name)
    actives, inactives = ds.get_classifications_overview()
    data  = {
        'name': dataset_name,
        'actives': int(actives),
        'inactives': int(inactives),
        'set_type': ds.set_type,
        'compounds': ds.compounds,
        'activities': ds.activities
    }
    return json.dumps(data)

@login_required
@app.route('/get_bioprofile/<profile_name>')
def get_bioprofile(profile_name):
    json_filename = os.path.join(g.user.get_user_folder('profiles'), '{}.json'.format(profile_name))
    with open(json_filename) as json_file:
        json_data = json.load(json_file)
    return json.dumps(json_data)


@login_required
@app.route('/download_bioprofile/<profile_name>')
def download_bioprofile(profile_name):
    bioprofile = g.user.load_bioprofile(profile_name)
    as_frame = bioprofile.to_frame()
    return as_frame.to_csv()

@login_required
@app.route('/get_bioprofile_class_overview/<profile_name>')
def get_bioprofile_class_overview(profile_name):
    bioprofile = g.user.load_bioprofile(profile_name)

    class_overview = bioprofile.classification_overview().to_dict()

    # convert to a json format for plotting in d3
    json_data = []

    for aid in class_overview['inactives']:
        row = {}
        row['inactives'] = int(class_overview['inactives'][aid])
        row['actives'] = int(class_overview['actives'][aid])
        row['aid'] = int(aid)
        json_data.append(row)

    return json.dumps(json_data)

@login_required
@app.route('/get_bioprofile_descriptions/<profile_name>')
def get_bioprofile_descriptions(profile_name):

    bioprofile = g.user.load_bioprofile(profile_name)

    class_overview = bioprofile.classification_overview()
    descriptions = bioprofile.get_bioassay_info()


    for desc in descriptions:
        desc['no_inactives'] = int(class_overview.loc[int(desc['AID']), 'inactives'])
        desc['no_actives'] = int(class_overview.loc[int(desc['AID']), 'actives'])

        # the desciption is too long so just take the first sentence
        desc['Description'] = desc['Description']
        desc['AID'] = int(desc['AID'])


    return json.dumps(descriptions)



@login_required
@app.route('/delete_profile', methods=['POST'])
def delete_profile():
    """
        a remove a specific bioprofile

    """
    json_data = request.get_json()

    profile_name = json_data['profile_name']

    profile_path = os.path.join(g.user.get_user_folder('profiles'), '{}.json'.format(profile_name))

    os.remove(profile_path)

    return 'OK', 200

@login_required
@app.route('/delete_dataset', methods=['POST'])
def delete_dataset():
    """
        a remove a specific bioprofile

    """
    json_data = request.get_json()

    dataset_name = json_data['dataset_name']

    dataset_path = os.path.join(g.user.get_user_folder('datasets'), '{}.json'.format(dataset_name))

    os.remove(dataset_path)

    return 'OK', 200



@login_required
@app.route('/sendClusterData', methods=['POST'])
def send_cluster_data():
    """
        retrieve cluster data from the server

    """
    json_data = request.get_json()

    cluster_data = json_data['results']

    fp_profile = g.user.load_fp_profile(cluster_data['currentClustering'])

    profile = g.user.load_bioprofile(fp_profile.meta['profile_used'])

    clusters = pd.DataFrame(cluster_data['clusterAssignments'])

    grouped_by_cluster = clusters.groupby('classLabel')

    main_frame = pd.DataFrame(list(zip(profile.cids, profile.outcomes)), index=profile.aids, columns=['cids', 'outcomes'])
    main_frame.index = main_frame.index.astype('int')
    main_frame.cids = main_frame.cids.astype('int')
    main_frame.outcomes = main_frame.outcomes.astype('int')

    stats_frame = pd.DataFrame(profile.stats)
    stats_frame.index = stats_frame.aid.astype(int)


    conversions = {"Specificity": "float", "Coverage": "float", "FP": "float", "CCR": "float",
                   "L parameter": "float", "Sensitivity": "float", "PPV": "float", "NPV": "float",
                   "TP": "int", "FN": "int", "aid": "int", "TN": "int"}

    for col, conversion in conversions.items():
        stats_frame[col] = stats_frame[col].astype(eval(conversion))

    # remove previous clusters

    previous_files = glob.glob(os.path.join(g.user.get_user_folder('profiles'), '{}_cluster_*.json'.format(profile.name)))

    for json_file in previous_files:
        os.remove(json_file)

    sorted_groups = sorted(grouped_by_cluster, key=lambda data: len(data[1]), reverse=True)

    new_cluster_id = 0

    for clst, clstr_data in sorted_groups:
        clstr_aids = clstr_data.aid.tolist()

        sub_frame = main_frame.loc[clstr_aids]
        sub_stats = stats_frame.loc[clstr_aids]

        for col, conversion in conversions.items():
            sub_stats[col] = sub_stats[col].astype(eval(conversion))


        aids = sub_frame.index.astype(int).tolist()
        cids = sub_frame.cids.astype(int).tolist()
        outcomes = sub_frame.outcomes.astype(int).tolist()
        stats = sub_stats.to_dict('records')
        name = '{}_cluster_{}'.format(profile.name, new_cluster_id)

        meta = {"num_cmps": len(cids),
                "num_total_actives": int((sub_frame['outcomes'] == 1).sum()),
                "num_total_inactives": int((sub_frame['outcomes'] == -1).sum()),
                "training_set": profile.meta['training_set'],
                "num_aids": len(aids)}

        sub_profile = bp.Bioprofile(name, cids, aids, outcomes, stats, meta)

        sub_profile.to_json(g.user.get_user_folder('profiles'))

        new_cluster_id = new_cluster_id + 1

    return 'OK', 200


@login_required
@app.route('/get_adj_matrix/<clustering_name>')
def get_adj_matrix(clustering_name):

    adj_matrix = g.user.load_adj_matrix(clustering_name)

    data  = {
        'links': adj_matrix.links,
        'nodes': adj_matrix.nodes,
        'profile_used': adj_matrix.profile_used,
        'linkage': adj_matrix.linkage
    }


    return json.dumps(data)


@login_required
@app.route('/add_inhouse_dataset', methods=['POST'])
def add_inhouse_dataset():
    json_data = request.get_json()

    inhouse_ds_name = json_data['inhouse_ds_name']

    compounds, input_type = inhouse_databases.get_inhouse_database(inhouse_ds_name)

    identifiers = [compound['identifier'] for compound in compounds]
    activities = [compound['activity'] for compound in compounds]

    name = inhouse_ds_name.lower().replace(' ', '_')

    ds_io.write_ds_to_json(identifiers, activities, g.user.get_user_folder('datasets'), name, input_type, set_type='training')

    return 'OK', 200


@login_required
@app.route('/filter_profile', methods=['POST'])
def filter_profile():
    """
        a post method to receive a users profile filter
        choices in the CIIPro Optimizer tab

    """
    json_data = request.get_json()

    bioprofile = g.user.load_bioprofile(json_data['profile_name'])

    assays_to_drop = []
    for i, profile_stat in enumerate(bioprofile.stats):
        for filter in json_data['filters']:
            if profile_stat[filter['stat']] < float(filter['thresh']):
                assays_to_drop.append(profile_stat['aid'])

    # get all indices for each aid

    assays_to_drop = list(set(assays_to_drop))

    idxs = []
    for aid in set(assays_to_drop):
        idxs = idxs + [i for i, x in enumerate(bioprofile.aids) if x == aid]

    # now remove them from cids, aids and outcomes

    for idx in sorted(idxs, reverse=True):
        idxs
        del bioprofile.aids[idx]
        del bioprofile.cids[idx]
        del bioprofile.outcomes[idx]

    bioprofile.stats = [profile_stat for profile_stat in bioprofile.stats
                        if profile_stat['aid'] not in assays_to_drop]

    profile_matrix = bioprofile.to_frame()

    bioprofile.meta['num_total_actives'] = int((profile_matrix == 1).sum().sum())
    bioprofile.meta['num_total_inactives'] = int((profile_matrix == -1).sum().sum())
    bioprofile.meta['num_cmps'] = int(profile_matrix.shape[0])
    bioprofile.meta['num_aids'] = int(profile_matrix.shape[1])


    bioprofile.to_json(g.user.get_user_folder('profiles'))

    return 'OK', 200



if __name__ == '__main__': #says if this scripts is run directly, start the application
	app.run()
