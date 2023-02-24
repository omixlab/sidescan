from sidescan.cli import train_model_argument_parser
from sklearn.model_selection import train_test_split
from sklearn.multioutput import MultiOutputClassifier
from sklearn.metrics import classification_report, accuracy_score, precision_score, recall_score, roc_auc_score, f1_score
from sklearn.ensemble import ExtraTreesClassifier
from xgboost import XGBClassifier
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import TomekLinks, RandomUnderSampler
from flaml import AutoML
from tqdm import tqdm
import warnings
import pandas as pd
import numpy as np
import pickle
import json
import os

warnings.simplefilter('ignore')

def main():

    arguments = train_model_argument_parser.parse_args()

    df_features = pd.read_csv(os.path.join(arguments.directory, 'fingerprints.csv'))
    df_labels   = pd.read_csv(os.path.join(arguments.directory, 'labels.csv'))

    models_results = []

    models = {}

    f1_scores = []

    for column in tqdm(df_labels.columns):
        
        y = df_labels[column]

        X_train, X_test, y_train, y_test = train_test_split(df_features, y)

        try:
            X_train, y_train = TomekLinks().fit_resample(X_train, y_train)
            X_train, y_train = SMOTE().fit_resample(X_train, y_train)
            X_test, y_test = RandomUnderSampler().fit_resample(X_test, y_test)
        except:
            print('unable to train model for', column)
        models[column] = ExtraTreesClassifier() #AutoML(task='classification', estimator_list=["lgbm"], time_budget=60, verbose=0)
        models[column].fit(X_train, y_train)

        y_pred = models[column].predict(X_test)

        try:
            f1 = f1_score(y_test, y_pred)
        except:
            f1 = 0

        try:
            roc_auc = roc_auc_score(y_test, y_pred)
        except:
            roc_auc = 0

        try:
            recall = recall_score(y_test, y_pred)
        except:
            recall = 0

        try:
            precision = precision_score(y_test, y_pred)
        except:
            precision = 0

        try:
            accuracy = accuracy_score(y_test, y_pred)
        except:
            accuracy = 0

        models[column].f1 = f1
        models[column].roc_auc = roc_auc
        models[column].recall = recall
        models[column].precision = precision
        models[column].accuracy = accuracy

        models_results.append({'effect': column, 'f1': f1, 'roc_auc': roc_auc, 'recall': recall, 'precision': precision, 'accuracy': accuracy})

    with open(os.path.join(arguments.directory, 'models.pickle'), 'wb') as writer:
        writer.write(pickle.dumps(models))

    df_models_results = pd.DataFrame(models_results).sort_values('f1', ascending=False)
    df_models_results.to_csv(os.path.join(arguments.directory, 'models.scores.csv'))

if __name__ == '__main__':
    main()
    
