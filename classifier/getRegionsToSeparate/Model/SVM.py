# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from sklearn import svm
from joblib import dump

def SVM(X_matrix, y_array):
    X = X_matrix
    y = y_array
    clf = svm.SVC(gamma='scale')
    clf.fit(X, y)
    dump(clf, './Model/models/svm.joblib')
