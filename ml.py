from sklearn.metrics import f1_score, accuracy_score
from sklearn.metrics import cohen_kappa_score, matthews_corrcoef, precision_score, recall_score, confusion_matrix


def get_class_stats(predicted_classes, y):


    acc = accuracy_score(y, predicted_classes)
    f1_sc = f1_score(y, predicted_classes)
    cohen_kappa = cohen_kappa_score(y, predicted_classes)
    matthews_corr = matthews_corrcoef(y, predicted_classes)
    precision = precision_score(y, predicted_classes)
    recall = recall_score(y, predicted_classes)

    # Specificity calculation
    tn, fp, fn, tp = confusion_matrix(y, predicted_classes).ravel()
    specificity = tn / (tn + fp)

    return {'ACC': float(acc), 'F1Score': float(f1_sc), 'Cohens Kappa': float(cohen_kappa),
            'MCC': float(matthews_corr), 'Precision': float(precision),
            'Recall': float(recall), 'Specificity': float(specificity)}
