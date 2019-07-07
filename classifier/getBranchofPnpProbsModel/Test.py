#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 13:15:34 2019

@author: kuangmeng
"""
from joblib import load
import sys
methods = ["randomforest", "svm", "abc"]
from data_generate import Prepare_Train_Data

def Generate_Data(file1, file2):
    train_data, train_label = Prepare_Train_Data(file1, file2)
    return train_data, train_label

def Test_L(path, test_data):
    clf = load(path)
    result = clf.predict(test_data)
    return int(result[0])

def Test(test_data, flag):
    l_path = "./Models/models/"
    methods = ["randomforest", "svm", "sgd", "abc"]
    tail = ".joblib"
    path = l_path + methods[int(flag)] + tail
    result = Test_L(path, test_data)
    return result


def getAccuracy(test, ref):
    lens = len(test)
    tmp_num= 0
    for i in range(lens):
        if test[i] == ref[i]:
            tmp_num += 1
    return tmp_num / float(lens) 


def TestFile(train, label):
    test = train[0: int(len(train) * 0.3)]
    test_label = label[0: int(len(train) * 0.3)]
    total_result = []
    for flag in range(3):
        result = []
        for test_i in test:
            result.append(Test([test_i], flag))
        total_result.append(result)
    
    for i in range(3):
        print(methods[i])
        print(getAccuracy(total_result[i], test_label))

if __name__ == "__main__":
    train, label = Generate_Data("./data/pnp1.txt", "./data/pnp2.txt")
    TestFile(train, label)

