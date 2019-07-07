#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 20:56:25 2019

@author: kuangmeng
"""

import numpy as np

def getPara(para, file):
    filein = open(file, "r").read().splitlines()
    max_avg = 0
    min_avg = 0
    max_sd = 0
    min_sd = 0
    max_blo = 0
    min_blo = 0
    if len(para) != 0 :
        max_avg = para[0]
        min_avg = para[1]
        max_sd = para[2]
        min_sd = para[3]
        max_blo = para[4]
        min_blo = para[5]
    for i in range(len(filein)):
        tmp_str = filein[i].split("\t")
        if float(tmp_str[0]) > max_avg:
            max_avg = float(tmp_str[0])
        elif float(tmp_str[0]) < min_avg:
            min_avg = float(tmp_str[0])
        if float(tmp_str[1]) > max_sd:
            max_sd = float(tmp_str[1])
        elif float(tmp_str[1]) < min_sd:
            min_sd = float(tmp_str[1])
        if float(tmp_str[2]) > max_blo:
            max_blo = float(tmp_str[2])
        elif float(tmp_str[2]) < min_blo:
            min_blo = float(tmp_str[2])
    return [max_avg, min_avg, max_sd, min_sd, max_blo, min_blo]

def getTrain(file, flag, para):
    traindata = []
    trainlabel = []
    num = 0
    filein = open(file, "r").read().splitlines()
    line_num = int(len(filein) * 0.7)
    max_avg = 0
    min_avg = 0
    max_sd = 0
    min_sd = 0
    max_blo = 0
    min_blo = 0
    if len(para) != 0 :
        max_avg = para[0]
        min_avg = para[1]
        max_sd = para[2]
        min_sd = para[3]
        max_blo = para[4]
        min_blo = para[5]
    else:
        for i in range(len(filein)):
            tmp_str = filein[i].split("\t")
            if float(tmp_str[0]) > max_avg:
                max_avg = float(tmp_str[0])
            elif float(tmp_str[0]) < min_avg:
                min_avg = float(tmp_str[0])
            if float(tmp_str[1]) > max_sd:
                max_sd = float(tmp_str[1])
            elif float(tmp_str[1]) < min_sd:
                min_sd = float(tmp_str[1])
            if float(tmp_str[2]) > max_blo:
                max_blo = float(tmp_str[2])
            elif float(tmp_str[2]) < min_blo:
                min_blo = float(tmp_str[2])

    for i in range(line_num):
        tmp_str = filein[i].split("\t")
        if float(tmp_str[0]) == 0 and float(tmp_str[1]) == 0 and float(tmp_str[2]) == 0:
            continue
        else:
            num += 1
            traindata.append([(float(tmp_str[0]) - min_avg)/ (max_avg - min_avg), (float(tmp_str[1]) - min_sd) / (max_sd - min_sd), (float(tmp_str[2]) - min_blo) / (max_blo - min_blo)])
            trainlabel.append([int(flag)])

    return traindata, trainlabel

def gen_data(file, data, label, count):
    with open(file,"w") as f:
        f.write("%d, 3\n" % count)
        #产生一个随机坐标(x,y)
        for i in range(count):
            x = float(data[i][0])
            y = float(data[i][1])
            z = float(data[i][2])
            if len(label) == 0:
                t = 0
            else:
                t = int(label[i])
            f.write("%f, %f, %f, %d\n" % (x, y, z, t))

def Prepare_Train_Data(file1, file2):
    para = getPara([], file1)
    para = getPara(para, file2)
    train_data1, train_label1 = getTrain(file1, 0, para)
    train_data2, train_label2 = getTrain(file2, 1, para)
    train_data = train_data1 + train_data2
    train_label = np.append(train_label1, train_label2)
    gen_data("./tmp/train.csv", train_data, train_label, len(train_label))
    return train_data, train_label
