#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:18:14 2019

@author: mmkuang
"""
import math
import time

tmp_str = "ARNDCQEGHILKMFPSTWYV"

matrix = []

def getMatrix():
    matrix.append([4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0])
    matrix.append([-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3])
    matrix.append([-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3])
    matrix.append([-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3])
    matrix.append([0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1])
    matrix.append([-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2])
    matrix.append([-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2])
    matrix.append([0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3])
    matrix.append([-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3])
    matrix.append([-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3])
    matrix.append([-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1])
    matrix.append([-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2])
    matrix.append([-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1])
    matrix.append([-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1])
    matrix.append([-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2])
    matrix.append([1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2])
    matrix.append([0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0])
    matrix.append([-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3])
    matrix.append([-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1])
    matrix.append([0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4])

def calculateColScore(real_pnp):
    getMatrix()
    col_score = []
    peak_length_ratio = 0.0
    file_context = real_pnp.split("\n")
    dic = {}
    has_key = False
    value = ""
    key = ""
    for itm in range(len(file_context)):
        if file_context[itm][0:1] == ">":
            if has_key == True:
                dic[key] = value
                value = ""
                key = ""
                has_key = False
            has_key = True
            key = file_context[itm]
        elif has_key == True:
            value = value.replace("\r","") + file_context[itm].replace("\r","")
    dic[key] = value
    dickeys =sorted(dic.keys())
    lens_ = (len(dickeys) * (len(dickeys) - 1)) / 2
    lens = len(value)
    tmp_un_sp = 0.0
    for i in range(lens):
        tmp_score = 0.0
        for k1 in range(len(dickeys) - 1):
            for k2 in range(k1 + 1, len(dickeys)):
                if getIdx(dic[dickeys[k1]][i]) == -1 or getIdx(dic[dickeys[k2]][i]) == -1:
                    tmp_score += 0
                else:
                    tmp_score += float(matrix[getIdx(dic[dickeys[k1]][i])][getIdx(dic[dickeys[k2]][i])])
        tmp_score /= lens_
        tmp_un_sp += tmp_score
        col_score.append(tmp_score)
    print("[SUPPORT STEP] Calculated Column Score!")
    tmp_sd = 0.0
    if lens != 0:
        tmp_un_sp /= lens
        tmp_sd = getSD(tmp_un_sp, col_score, lens)
        peak_length_ratio = getPeakLengthRatio(col_score, lens)
    else:
        tmp_un_sp = 0
    prepare_data_2 = time.time()
    return prepare_data_2, col_score, tmp_un_sp, lens, len(dickeys), tmp_sd, peak_length_ratio

def getAvgColScore(real_msa):
    file_context = []
    with open(real_msa, 'r') as filein:
        file_context = filein.read().splitlines()
    getMatrix()
    dic = {}
    has_key = False
    value = ""
    key = ""
    for itm in range(len(file_context)):
        if file_context[itm][0:1] == ">":
            if has_key == True:
                dic[key] = value
                value = ""
                key = ""
                has_key = False
            has_key = True
            key = file_context[itm]
        elif has_key == True:
            value = value.replace("\r","") + file_context[itm].replace("\r","")
    dic[key] = value
    dickeys =sorted(dic.keys())
    lens_ = (len(dickeys) * (len(dickeys) - 1)) / 2
    lens = len(value)
    if lens_ * lens == 0:
        return -1
    tmp_un_sp = 0
    for i in range(lens):
        tmp_score = 0.0
        for k1 in range(len(dickeys) - 1):
            for k2 in range(k1 + 1, len(dickeys)):
                if getIdx(dic[dickeys[k1]][i]) == -1 or getIdx(dic[dickeys[k2]][i]) == -1:
                    tmp_score += 0
                else:
                    tmp_score += float(matrix[getIdx(dic[dickeys[k1]][i])][getIdx(dic[dickeys[k2]][i])])
        tmp_score /= lens_
        tmp_un_sp += tmp_score
    return float(tmp_un_sp / lens)

def getSD(tmp_un_sp, tmp_un_sp_arr, lens):
    sd_ = 0.0
    for item in tmp_un_sp_arr:
        sd_ += (float(item) - float(tmp_un_sp)) ** 2
    sd_ /= lens
    return math.sqrt(sd_)

def getPeakLengthRatio(tmp_un_sp_arr, lens):
    ratio = 0.0
    for item in tmp_un_sp_arr:
        if float(item) >= 1.0:
            ratio += 1
    return ratio / lens

def getIdx(char):
    for i in range(len(tmp_str)):
        if char == tmp_str[i]:
            return i
    return -1
