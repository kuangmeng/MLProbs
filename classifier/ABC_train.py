#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 13:11:36 2018
@author: kuangmeng
"""

from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier
from joblib import dump

def ABC(X_matrix, y_array):
    X = X_matrix
    Y = y_array
    model = DecisionTreeClassifier(max_depth=5)
    clf = AdaBoostClassifier(model, n_estimators=200)
    clf.fit(X, Y)
    dump(clf, './model/abc.joblib')
