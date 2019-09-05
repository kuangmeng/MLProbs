#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 13:15:34 2019

@author: kuangmeng
"""
from joblib import load
import sys
from script import PrepareData

methods = ["randomforest", "svm", "abc"]

def Generate_Data():
    train_data, train_label = PrepareData()
    return train_data, train_label

def Test_L(path, test_data):
    clf = load(path)
    result = clf.predict(test_data)
    return int(result[0])

def Test(test_data, flag):
    l_path = "./Model/models/"
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
    test = train[int(len(train)*0.7): ]
    test_label = label[int(len(train)*0.7): ]
    total_result = []
    for flag in range(3):
        result = []
        for test_i in test:
            result.append(Test([test_i], flag))
        total_result.append(result)
    
    for i in range(3):
    #     print(methods[i])
        print(getAccuracy(total_result[i], test_label))
    #     print(total_result[i])


if __name__ == "__main__":
    train, label = Generate_Data()
    TestFile(train, label)
