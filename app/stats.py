from sklearn.metrics import cohen_kappa_score, matthews_corrcoef, precision_score, recall_score, confusion_matrix
from sklearn.metrics import make_scorer
from sklearn.metrics import roc_auc_score, roc_curve, auc, f1_score, accuracy_score, r2_score, \
    max_error, \
    mean_squared_error, mean_absolute_percentage_error, d2_pinball_score


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
            'MCC': matthews_corr, 'Precision/PPV': precision, 'Recall': recall, 'Specificity': specificity, 'CCR': ccr,
            'R2-score': None, 'Explained-variance': None, 'Max-error': None, 'Mean-squared-error': None,
            'Mean-absolute-percentage-error': None, 'D2-pinball-score': None
            }


# scoring dictionary, just a dictionary containing the evaluation metrics passed through a make_scorer()
# fx, necessary for use in GridSearchCV

class_scoring = {'ACC': make_scorer(accuracy_score), 'F1-Score': make_scorer(f1_score),
                 'AUC': make_scorer(roc_auc_score),
                 'Cohen\'s Kappa': make_scorer(cohen_kappa_score), 'MCC': make_scorer(matthews_corrcoef),
                 'Precision': make_scorer(precision_score), 'Recall': make_scorer(recall_score),
                 'R2-score': None, 'Max-error': None, 'Mean-squared-error': None,
                 'Mean-absolute-percentage-error': None, 'D2-pinball-score': None
                 }


def get_regress_stats(model, X, y):
    """

    :param model: If None, assume X == y_true and y == y_pred, else should be a trained model
    :param X: Data to predict
    :param y: correct value
    :return:
    """
    if not model:
        predicted_values = y.copy()
        y = X
    else:
        if 'predict_values' in dir(model):
            predicted_values = model.predict_values(X, verbose=0)[:, 0]
        else:
            predicted_values = model.predict(X)

    r2 = r2_score(y, predicted_values)
    m_error = max_error(y, predicted_values)
    ms_error = mean_squared_error(y, predicted_values)
    mape = mean_absolute_percentage_error(y, predicted_values)
    pinball_score = d2_pinball_score(y, predicted_values)

    return {'ACC': None, 'F1-Score': None, 'AUC': None, 'Cohen\'s Kappa': None,
            'MCC': None, 'Precision/PPV': None, 'Recall': None, 'Specificity': None, 'CCR': None,
            'R2-score': r2, 'Max-error': m_error, 'Mean-squared-error': ms_error,
            'Mean-absolute-percentage-error': mape, 'D2-pinball-score': pinball_score}


# scoring dictionary, just a dictionary containing the evaluation metrics passed through a make_scorer()
# fx, necessary for use in GridSearchCV

regress_scoring = {'Explained-variance': make_scorer(r2_score), 'Max-error': make_scorer(max_error),
                   'Mean-squared-error': make_scorer(mean_squared_error),
                   'Mean-absolute-percentage-error': make_scorer(mean_absolute_percentage_error),
                   'D2-pinball-score': make_scorer(d2_pinball_score)}
