from sklearn.metrics import cohen_kappa_score, matthews_corrcoef, precision_score, recall_score, confusion_matrix
from sklearn.metrics import make_scorer
from sklearn.metrics import roc_auc_score, roc_curve, auc, f1_score, accuracy_score
from sklearn.metrics import r2_score, mean_absolute_error, mean_absolute_percentage_error, \
    mean_squared_error, mean_squared_log_error, explained_variance_score, median_absolute_error, mean_tweedie_deviance


def regression_roc_auc_score(y_true, y_pred, num_rounds=10000):
    """
    Computes Regression-ROC-AUC-score.

    Parameters:
    ----------
    y_true: array-like of shape (n_samples,). Binary or continuous target variable.
    y_pred: array-like of shape (n_samples,). Target scores.
    num_rounds: int or string. If integer, number of random pairs of observations.
                If string, 'exact', all possible pairs of observations will be evaluated.

    Returns:
    -------
    rroc: float. Regression-ROC-AUC-score.
    """

    import numpy as np

    y_true = np.array(y_true)
    y_pred = np.array(y_pred)

    num_pairs = 0
    num_same_sign = 0

    for i, j in _yield_pairs(y_true, num_rounds):
        diff_true = y_true[i] - y_true[j]
        diff_score = y_pred[i] - y_pred[j]
        if diff_true * diff_score > 0:
            num_same_sign += 1
        elif diff_score == 0:
            num_same_sign += .5
        num_pairs += 1

    return num_same_sign / num_pairs


def _yield_pairs(y_true, num_rounds):
    """
    Returns pairs of valid indices. Indices must belong to observations having different values.

    Parameters:
    ----------
    y_true: array-like of shape (n_samples,). Binary or continuous target variable.
    num_rounds: int or string. If integer, number of random pairs of observations to return.
                If string, 'exact', all possible pairs of observations will be returned.

    Yields:
    -------
    i, j: tuple of int of shape (2,). Indices referred to a pair of samples.

    """
    import numpy as np

    if num_rounds == 'exact':
        for i in range(len(y_true)):
            for j in np.where((y_true != y_true[i]) & (np.arange(len(y_true)) > i))[0]:
                yield i, j
    else:
        for r in range(num_rounds):
            i = np.random.choice(range(len(y_true)))
            j = np.random.choice(np.where(y_true != y_true[i])[0])
            yield i, j


def get_class_stats(model, X, y):
    """

    :param model: If None, assume X == y_true and y == y_pred, else should be a trained model
    :param X: Data to predict
    :param y: correct classes
    :return:
    """
    if not model:
        predicted_probas = y
        predicted_classes = y.copy()
        predicted_classes[predicted_classes >= 0.5] = 1
        predicted_classes[predicted_classes < 0.5] = 0
        y = X
    else:
        if 'predict_classes' in dir(model):
            predicted_classes = model.predict_classes(X, verbose=0)[:, 0]
            predicted_probas = model.predict_proba(X, verbose=0)[:, 0]
        else:
            predicted_classes = model.predict(X)
            predicted_probas = model.predict_proba(X)[:, 1]

    acc = accuracy_score(y, predicted_classes)
    f1_sc = f1_score(y, predicted_classes)

    # Sometimes SVM spits out probabilties with of inf
    # so set them as 1
    from numpy import inf
    predicted_probas[predicted_probas == inf] = 1

    fpr_tr, tpr_tr, thresholds_tr = roc_curve(y, predicted_probas)
    roc_auc = auc(fpr_tr, tpr_tr)
    # test classification results

    cohen_kappa = cohen_kappa_score(y, predicted_classes)
    matthews_corr = matthews_corrcoef(y, predicted_classes)
    precision = precision_score(y, predicted_classes)
    recall = recall_score(y, predicted_classes)

    # Specificity calculation
    tn, fp, fn, tp = confusion_matrix(y, predicted_classes).ravel()
    specificity = tn / (tn + fp)
    ccr = (recall + specificity) / 2

    return {'ACC': acc, 'F1-Score': f1_sc, 'AUC': roc_auc, 'Cohen\'s Kappa': cohen_kappa,
            'MCC': matthews_corr, 'Precision/PPV': precision, 'Recall': recall, 'Specificity': specificity, 'CCR': ccr}


def get_regressive_stats(model, X, y):
    """

    :param model: If None, assume X == y_true and y == y_pred, else should be a trained model
    :param X: Data to predict
    :param y: correct classes
    :return:
    """
    if not model:
        predicted_probas = y
        y = X
    else:
        if 'predict_classes' in dir(model):
            predicted_values = model.predict_classes(X, verbose=0)[:, 0]
        else:
            predicted_values = model.predict(X)

    r2 = r2_score(y, predicted_values)
    mae = mean_absolute_error(y, predicted_values)
    mape = mean_absolute_percentage_error(y, predicted_values)
    mse = mean_squared_error(y, predicted_values)
    msle = mean_squared_log_error(y, predicted_values)
    medae = median_absolute_error(y, predicted_values)

    roc_auc = regression_roc_auc_score(y, predicted_values)

    return {'R2': r2, 'mean-absolute-error': mae, 'AUC': roc_auc, 'mean-absolute-percentage-error': mape,
            'mean-square-error': mse,'mean-square-log-error': msle, 'median-absolute-error': medae}


# scoring dictionary, just a dictionary containing the evaluation metrics passed through a make_scorer()
# fx, necessary for use in GridSearchCV

class_scoring = {'ACC': make_scorer(accuracy_score), 'F1-Score': make_scorer(f1_score),
                 'AUC': make_scorer(roc_auc_score),
                 'Cohen\'s Kappa': make_scorer(cohen_kappa_score), 'MCC': make_scorer(matthews_corrcoef),
                 'Precision': make_scorer(precision_score), 'Recall': make_scorer(recall_score)}

regressive_scoring = {'R2': make_scorer(r2_score), 'Mean_absolute_error': make_scorer(mean_absolute_error),
                      'AUC': make_scorer(regression_roc_auc_score),
                      'Mean_absolute_percentage_error': make_scorer(mean_absolute_percentage_error),
                      'Mean_square_error': make_scorer(mean_squared_error),
                      'Mean_square_log_error': make_scorer(mean_squared_log_error),
                      'Median_absolute_error': make_scorer(median_absolute_error)}
