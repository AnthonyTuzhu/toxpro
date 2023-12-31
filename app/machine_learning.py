import numpy as np
import pandas as pd

from sklearn import model_selection
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, AdaBoostClassifier, AdaBoostRegressor
from sklearn.naive_bayes import GaussianNB, BernoulliNB
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.svm import SVC, SVR
from sklearn import model_selection
from sklearn import pipeline
from sklearn.model_selection import cross_val_predict
from sklearn.preprocessing import StandardScaler

import pickle

from app.stats import class_scoring, regress_scoring, get_class_stats, get_regress_stats

seed = 0

CLASSIFIER_ALGORITHMS = [
    ('RF', RandomForestClassifier(max_depth=10,  # max depth 10 to prevent overfitting
                                  class_weight='balanced',
                                  random_state=seed), {'RF__n_estimators': [5, 10, 25, 100, 200]}),
    ('kNN', KNeighborsClassifier(metric='euclidean'), {'kNN__n_neighbors': [1, 3, 5, 10, 20],
                                                       'kNN__weights': ['uniform', 'distance']}),
    ('SVM', SVC(probability=True,
                class_weight='balanced',
                random_state=seed), {'SVM__kernel': ['linear', 'rbf', 'poly'],
                                     'SVM__gamma': [1e-2, 1e-3],
                                     'SVM__C': [0.1, 1, 10]}),
    ('BNB', BernoulliNB(alpha=1.0), {}),
    ('ADA', AdaBoostClassifier(n_estimators=100, learning_rate=0.9, random_state=seed), {})
]

CLASSIFIER_ALGORITHMS_DICT = dict([(t[0], (t[0], t[1], t[2])) for t in CLASSIFIER_ALGORITHMS])

REGRESSOR_ALGORITHMS = [
    ('RF', RandomForestRegressor(max_depth=10,  # max depth 10 to prevent overfitting
                                 random_state=seed), {'RF__n_estimators': [5, 10, 25, 100, 200]}),
    ('kNN', KNeighborsRegressor(metric='euclidean'), {'kNN__n_neighbors': [1, 3, 5, 10, 20],
                                                      'kNN__weights': ['uniform', 'distance']}),
    ('SVR', SVR(), {'SVM__kernel': ['linear', 'rbf', 'poly'],
                    'SVM__gamma': [1e-2, 1e-3],
                    'SVM__C': [0.1, 1, 10]}),
    ('ADA', AdaBoostRegressor(n_estimators=100, learning_rate=0.9, random_state=seed), {})
]

REGRESSOR_ALGORITHMS_DICT = dict([(t[0], (t[0], t[1], t[2])) for t in REGRESSOR_ALGORITHMS])


def build_qsar_model(X: pd.DataFrame, y: pd.Series, alg: str, scale=True):
    """ build a classification QSAR model """
    cv = model_selection.StratifiedKFold(shuffle=True, n_splits=5, random_state=seed)
    model_entry = CLASSIFIER_ALGORITHMS_DICT[alg]

    name, model, params = model_entry

    if scale:
        pipe = pipeline.Pipeline([('scaler', StandardScaler()), (name, model)])
    else:
        pipe = pipeline.Pipeline([(name, model)])

    grid_search = model_selection.GridSearchCV(pipe,
                                               param_grid=params,
                                               cv=cv,
                                               scoring=class_scoring,
                                               refit='AUC')
    grid_search.fit(X, y)
    best_estimator = grid_search.best_estimator_

    # get the predictions from the best performing model in 5 fold cv
    cv_predictions = pd.DataFrame(
        cross_val_predict(best_estimator, X, y,
                          cv=cv,
                          method='predict_proba'),
        index=y.index)

    cv_class = cv_predictions[1].copy()
    cv_class[cv_class >= 0.5] = 1
    cv_class[cv_class < 0.5] = 0
    five_fold_stats = get_class_stats(None, y, cv_predictions[1])

    # record the predictions and the results
    final_cv_predictions = pd.concat([cv_predictions[1], cv_class], axis=1)

    binary_model = pickle.dumps(best_estimator)

    # train_stats = get_class_stats(best_estimator, X, y)

    return binary_model, final_cv_predictions, five_fold_stats


def build_qsar_model_regression(X: pd.DataFrame, y: pd.Series, alg: str, scale=True):
    """ build a regressive QSAR model """
    cv = model_selection.StratifiedKFold(shuffle=True, n_splits=5, random_state=seed)
    model_entry = REGRESSOR_ALGORITHMS_DICT[alg]

    name, model, params = model_entry

    if scale:
        pipe = pipeline.Pipeline([('scaler', StandardScaler()), (name, model)])
    else:
        pipe = pipeline.Pipeline([(name, model)])

    grid_search = model_selection.GridSearchCV(pipe,
                                               param_grid=params,
                                               cv=cv,
                                               scoring=regress_scoring,
                                               refit='Mean-squared-error')
    grid_search.fit(X, y)
    best_estimator = grid_search.best_estimator_

    cv_predictions = pd.DataFrame(
        cross_val_predict(best_estimator, X, y,
                          cv=cv),
        index=y.index)

    five_fold_stats = get_regress_stats(None, y, cv_predictions)

    # record the predictions and the results
    final_cv_predictions = cv_predictions
    regressive_model = pickle.dumps(best_estimator)

    # train_stats = get_class_stats(best_estimator, X, y)

    return regressive_model, final_cv_predictions, five_fold_stats


if __name__ == '__main__':
    print(CLASSIFIER_ALGORITHMS_DICT)
