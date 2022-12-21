# this module is outlined
# here: https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xxii-background-jobs
from app import create_app
from app.db_models import db, Dataset, Chemical, QSARModel, Task
from app import chem_io
from app import machine_learning as ml

import sys, time
import datetime
from datetime import timezone
from rq import get_current_job
import pandas as pd

app = create_app()
app.app_context().push()

def build_qsar(user_id, dataset_name, descriptors, algorithm):

    job = get_current_job()
    job.meta['progress'] = 'Creating Features...'
    job.save_meta()
    query_statement = db.session.query(Chemical).join(Dataset,
                                                      Dataset.id == Chemical.dataset_id) \
        .filter(Dataset.dataset_name == dataset_name) \
        .filter(Dataset.user_id == user_id).statement
    df = pd.read_sql(query_statement, db.session.bind)

    dataset = Dataset.query.filter_by(dataset_name=dataset_name).first()
    # erase exists models
    qsar_model = QSARModel.query.filter_by(user_id=user_id,
                                           algorithm=algorithm,
                                           descriptors=descriptors,
                                           dataset_id=dataset.id).first()

    if qsar_model:
        db.session.delete(qsar_model)

    # create descriptors
    X = chem_io.get_desc(df, descriptors)

    y = df['activity']
    y.index = df['compound_id']
    y = y.loc[X.index]

    if descriptors == 'RDKit':
        scale = True
    else:
        scale = False

    # start training
    job.meta['progress'] = 'Training Model...'
    job.save_meta()
    model, cv_preds, train_stats = ml.build_qsar_model(X,
                                                       y,
                                                       algorithm,
                                                       scale=scale)


    qsar_model = QSARModel(user_id=user_id,
                           name=f'{dataset_name}-{descriptors}-{algorithm}',
                           algorithm=algorithm,
                           descriptors=descriptors,
                           dataset_id=dataset.id,
                           sklearn_model=model)



    # job done
    # think about putting this into its own helper function
    # like in _set_task_progress() here: https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xxii-background-jobs
    job.meta['progress'] = 'Complete'
    job.save_meta()
    task = Task.query.get(job.get_id())
    task.complete = True
    task.time_completed = datetime.datetime.now(timezone.utc)

    # not sure why I need to save this here
    # and not in the parent function
    db.session.add(qsar_model)
    db.session.commit()


def example(seconds):
    job = get_current_job()
    print('Starting task')
    for i in range(seconds):
        job.meta['progress'] = 100.0 * i / seconds
        job.save_meta()
        print(i)
        time.sleep(1)
    job.meta['progress'] = 100
    job.save_meta()
    print('Task completed')