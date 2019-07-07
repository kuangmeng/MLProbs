#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 13:18:30 2018

@author: kuangmeng
"""

from sklearn.ensemble import RandomForestClassifier
from joblib import dump

def RandomForest(X_matrix, y_array):
    X = X_matrix
    Y = y_array
    clf = RandomForestClassifier(n_estimators=10)
    clf = clf.fit(X, Y)
    dump(clf, './Models_Layer1/models/randomforest.joblib')
