# this module is outlined
# here: https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xxii-background-jobs
from app import create_app
from app.db_models import db, Dataset, Chemical, QSARModel, Task, CVResults
from app.curator.curator import Curator
from app import chem_io
from app import machine_learning as ml
from app import pubchem as pc

import sys, time
import datetime
from datetime import timezone
from rq import get_current_job
import pandas as pd

app = create_app()
app.app_context().push()


def build_qsar(user_id, dataset_name, descriptors, algorithm, type):
    """ function that ties together the entirity of building a QSAR model
     this is necessary to be able to submit as a task in redis. """

    job = get_current_job()
    # if Queue.started_job_registry is None:
    job.meta['progress'] = 'Creating Features...'
    job.save_meta()
    query_statement = db.session.query(Chemical).join(Dataset,
                                                      Dataset.id == Chemical.dataset_id) \
        .filter(Dataset.dataset_name == dataset_name) \
        .filter(Dataset.user_id == user_id).statement
    df = pd.read_sql(query_statement, db.session.connection())

    dataset = Dataset.query.filter_by(dataset_name=dataset_name, user_id=user_id).first()
    # erase exists models and results
    qsar_model = QSARModel.query.filter_by(user_id=user_id,
                                           algorithm=algorithm,
                                           descriptors=descriptors,
                                           type=type,
                                           dataset_id=dataset.id).first()

    if qsar_model:
        # erase exists models
        cv_results = CVResults.query.filter_by(qsar_model_id=qsar_model.id).first()

        if cv_results:
            db.session.delete(cv_results)
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
    if type == 'Classification':
        model, cv_preds, train_stats = ml.build_qsar_model(X,
                                                           y,
                                                           algorithm,
                                                           scale=scale)

    # Added regression models (share same results but parameters are different for each (using different columns)

    elif type == 'Regression':
        model, cv_preds, train_stats = ml.build_qsar_model_regression(X,
                                                                      y,
                                                                      algorithm,
                                                                      scale=scale)

    qsar_model = QSARModel(user_id=user_id,
                           name=f'{dataset_name}-{descriptors}-{algorithm}-{type}',
                           algorithm=algorithm,
                           descriptors=descriptors,
                           type=type,
                           dataset_id=dataset.id,
                           sklearn_model=model)
    db.session.add(qsar_model)
    db.session.commit()

    cv_results = CVResults(qsar_model_id=qsar_model.id,
                           accuracy=train_stats['ACC'],
                           f1_score=train_stats['F1-Score'],
                           area_under_roc=train_stats['AUC'],
                           cohens_kappa=train_stats['Cohen\'s Kappa'],
                           michaels_correlation=train_stats['MCC'],
                           precision=train_stats['Precision/PPV'],
                           recall=train_stats['Recall'],
                           specificity=train_stats['Specificity'],
                           correct_classification_rate=train_stats['CCR'],
                           r2_score=train_stats['R2-score'],
                           max_error=train_stats['Max-error'],
                           mean_squared_error=train_stats['Mean-squared-error'],
                           mean_absolute_percentage_error=train_stats['Mean-absolute-percentage-error'],
                           pinball_score=train_stats['D2-pinball-score']
                           )

    # job done think about putting this into its own helper function like in _set_task_progress() here:
    # https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xxii-background-jobs
    job.meta['progress'] = 'Complete'
    job.save_meta()
    task = Task.query.get(job.get_id())
    task.complete = True
    task.time_completed = datetime.datetime.now(timezone.utc)

    # not sure why I need to save this here
    # and not in the parent function

    db.session.add(cv_results)
    db.session.commit()


def curate_chems(user_id, dataset_name, duplicate_selection, replace=False):
    """ function that curates a set of chemicals submitted to redis """

    job = get_current_job()
    job.meta['progress'] = 'Curating...'
    job.save_meta()
    query_statement = db.session.query(Chemical).join(Dataset,
                                                      Dataset.id == Chemical.dataset_id) \
        .filter(Dataset.dataset_name == dataset_name) \
        .filter(Dataset.user_id == user_id).statement
    df = pd.read_sql(query_statement, db.session.connection())

    dataset = Dataset.query.filter_by(dataset_name=dataset_name, user_id=user_id).first()

    curator = Curator(df)
    curator.curate(duplicates=duplicate_selection)

    new_df = curator.new_df.copy()

    if not replace:

        job.meta['progress'] = 'Adding dataset...'
        job.save_meta()

        new_dataset = Dataset(dataset_name=dataset.dataset_name + '_curated', user_id=user_id)
        db.session.add(new_dataset)
        db.session.commit()

        for i, row in new_df.iterrows():
            activity = row['activity']
            cmp_id = row['compound_id']
            inchi = row['inchi']

            chem = Chemical(inchi=inchi, dataset_id=dataset.id, activity=activity, compound_id=cmp_id)
            new_dataset.chemicals.append(chem)

        db.session.add(new_dataset)
        db.session.commit()
    else:
        job.meta['progress'] = 'Replacing dataset...'
        job.save_meta()
        # remove all previous chemicals
        Chemical.query.filter_by(dataset_id=dataset.id).delete()

        db.session.commit()

        for i, row in new_df.iterrows():
            activity = row['activity']
            cmp_id = row['compound_id']
            inchi = row['inchi']

            chem = Chemical(inchi=inchi, dataset_id=dataset.id, activity=activity, compound_id=cmp_id)
            dataset.chemicals.append(chem)

        db.session.add(dataset)
        db.session.commit()

    # job done think about putting this into its own helper function like in _set_task_progress() here:
    # https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xxii-background-jobs
    job.meta['progress'] = 'Complete'
    job.save_meta()
    task = Task.query.get(job.get_id())
    task.complete = True
    task.time_completed = datetime.datetime.now(timezone.utc)

    # db.session.add(cv_results)
    # db.session.commit()


def add_pubchem_data(aid, user_id):
    """ function that adds data from pubchem """

    job = get_current_job()

    job.meta['progress'] = 'Importing data...'
    job.save_meta()
    df, fail_reason = pc.import_pubchem_aid(aid)

    if fail_reason is not None:
        job.meta['progress'] = f'Failed: {fail_reason}'
        job.save_meta()
        return

    job.meta['progress'] = 'Adding data...'
    job.save_meta()

    # check to see if user already imported a dataset
    # with this name
    Dataset.query.filter_by(user_id=user_id, dataset_name=f'AID_{aid}').delete()

    db.session.commit()

    dataset = Dataset(dataset_name=f'AID_{aid}', user_id=user_id)
    db.session.add(dataset)
    db.session.commit()

    for i, row in df.iterrows():
        activity = row['activity']
        cmp_id = row['compound_id']
        inchi = row['inchi']

        chem = Chemical(inchi=inchi, dataset_id=dataset.id, activity=activity, compound_id=cmp_id)
        dataset.chemicals.append(chem)

    db.session.add(dataset)
    db.session.commit()

    # job done
    # think about putting this into its own helper function
    # like in _set_task_progress() here: https://blog.miguelgrinberg.com/post/the-flask-mega-tutorial-part-xxii-background-jobs
    job.meta['progress'] = 'Complete'
    job.save_meta()
    task = Task.query.get(job.get_id())
    task.complete = True
    task.time_completed = datetime.datetime.now(timezone.utc)


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


if __name__ == '__main__':
    curate_chems(1, "AMES_curated")
