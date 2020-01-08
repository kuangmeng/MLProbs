#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:19:32 2019

@author: mmkuang
"""
from joblib import load

model_region = "./classifier/model/regions/randomforest.joblib"
para_region = "./classifier/model/regions/para.txt"

def getRealignStrategy(peak_length_ratio, avg_PID, sd_un_sp, un_sp):
    global para_region
    global model_region
    class_region = 1
    test = [peak_length_ratio, avg_PID, sd_un_sp, un_sp]
    para = []
    with open(para_region, 'r') as filein:
        file_context = filein.read().splitlines()
        for i in range(2 * len(test)):
            para.append(float(file_context[i]))
    real_test = []
    for i in range(len(test)):
        real_test.append((float(test[i]) - para[i * 2 + 1])/ (para[i * 2] - para[i * 2 + 1]))
    clf = load(model_region)
    class_region = clf.predict([real_test])[0]
    if class_region > 1 or class_region < 0:
        class_region = 1
    return class_region