#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:42:55 2019

@author: kuangmeng
"""

import sys
from script import PrepareData
sys.path.append("./Model/")

from SVM import SVM
from RandomForest import RandomForest
from ABC import ABC

def Generate_Data():
    train_data, train_label = PrepareData()
    return train_data, train_label

def Train(train_data, train_label):
    SVM(train_data, train_label)
    RandomForest(train_data, train_label)
    ABC(train_data, train_label)
    
if __name__ == "__main__":

    train_data, train_label = Generate_Data()
    Train(train_data, train_label)
