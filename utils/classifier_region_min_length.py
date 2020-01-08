#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:19:32 2019

@author: mmkuang
"""
from joblib import load
import time
model_lens = "./classifier/model/seq_lens/randomforest.joblib"
para_lens = "./classifier/model/seq_lens/para.txt"

def getRegionsLength(len_seqs, len_family, avg_PID, sd_PID, un_sp):
    global para_lens
    global model_lens
    class_lens = 3
    test = [len_seqs, len_family, avg_PID, sd_PID, un_sp]
    para = []
    with open(para_lens, 'r') as filein:
        file_context = filein.read().splitlines()
        for i in range(2 * len(test)):
            para.append(float(file_context[i]))
    real_test = []
    for i in range(len(test)):
        real_test.append((float(test[i]) - para[i * 2 + 1])/ (para[i * 2] - para[i * 2 + 1]))
    clf = load(model_lens)
    class_lens = clf.predict([real_test])[0]
    if class_lens > 3 or class_lens < 0:
        class_lens = 3
    return class_lens