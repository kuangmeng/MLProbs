#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:18:14 2019

@author: mmkuang
"""
tmp_str = "ARNDCQEGHILKMFPSTWYV"

matrix = [[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0] , [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3] , [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3] , [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3] , [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1] , [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2] , [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2] , [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3] , [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3] , [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3] , [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1] , [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2] , [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1] , [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1] , [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2] , [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2] , [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0] , [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3] , [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1] , [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],]

def detect_unreliable_regions(real_pnp, col_score):
    filein = open(real_pnp, 'r')
    file_context = filein.read().splitlines()
    filein.close()
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
    lens = len(value)
    fileout = open(col_score, 'w')
    fileout.write("#COL\t#COL_Score\n")
    for i in range(lens):
        tmp_score = 0.0
        for k1 in range(len(dickeys)):
            for k2 in range(k1, len(dickeys)):
                tmp_score += float(matrix[getIdx(dic[k1][i])][getIdx(dic[k2][i])])
        fileout.write(str(i) + "\t" + str(tmp_score) + "\n")
    fileout.write("#END\n")
    fileout.close()

def getIdx(char):
    for i in range(len(tmp_str)):
        if char == tmp_str[i]:
            return i
    return 0
