#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 16:54:37 2019

@author: mmkuang
"""
import os
from Detect_Unreliable_Family import Detect_Unreliable
def getUnreliableRegions(sigma, beta, theta, threshold, col_score, seq_file, real_output):
    # if Detect_Unreliable(theta, threshold, col_score, seq_file, real_output):
    #     return []
    last_col = 0
    num_score = []
    with open(col_score, 'r') as filein:
        file_context = filein.read().splitlines()
        for i in range(1, len(file_context) - 1):
            num_score.append(file_context[i].split())
            last_col = file_context[i].split()[0]
    unreliable_regions = []
    tmp_1 = 0
    tmp_2 = 0
    tmp_head = 0
    for item in num_score:
        if float(item[1]) <= sigma and float(item[1]) >= beta and tmp_1 == 0:
            tmp_head = int(item[0])
            tmp_1 = 1
        elif float(item[1]) <= sigma and float(item[1]) >= beta and tmp_1 == 1 and tmp_2 == 0:
            tmp_2 = 1
        elif float(item[1]) <= sigma and float(item[1]) >= beta and tmp_1 == 1 and tmp_2 == 1:
            if item[0] == last_col:
                #if (int(last_col) > 300 and int(item[0]) - int(tmp_head) > 15) or (int(last_col) < 300 and int(item[0]) - int(tmp_head) > 10):
                if int(item[0]) - int(tmp_head) > 15:
                    unreliable_regions.append([int(tmp_head), int(item[0]) - 1])
            else:
                continue
        elif (float(item[1]) > sigma or float(item[1]) < beta) and tmp_1 == 1 and tmp_2 == 1:
            #if (int(last_col) > 300 and int(item[0]) - int(tmp_head) > 15) or (int(last_col) < 300 and int(item[0]) - int(tmp_head) > 10):
            if int(item[0]) - int(tmp_head) > 15:
                unreliable_regions.append([int(tmp_head), int(item[0]) - 1])
            tmp_1 = 0
            tmp_2 = 0
            tmp_head = 0
        else:
            tmp_1 = 0
            tmp_2 = 0
            tmp_head = 0

    return unreliable_regions

def seperateUnreliableRegions(unreliable_regions, real_output, dir_output):
    filein = open(real_output, 'r')
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
    lens = len(value)
    dickeys = sorted(dic.keys())
    if len(unreliable_regions) == 0:
        os.system("cp " + real_output + " " + dir_output + "0-" + str(lens - 1) + ".reliable")
    else:
        if int(unreliable_regions[0][0]) > 1:
            fileout = open(dir_output + "0-" + str(int(unreliable_regions[0][0]) - 2) + ".reliable", 'w')
            for idx in range(len(dickeys)):
                fileout.write(dickeys[idx] + "\n")
                fileout.write(dic[dickeys[idx]][0: int(unreliable_regions[0][0]) - 1] + "\n")
            fileout.close()
        for item in unreliable_regions:
            fileout = open(dir_output + str(int(item[0]) - 1) + "-" + str(int(item[1]) - 1) + ".unreliable", 'w')
            for idx in range(len(dickeys)):
                fileout.write(dickeys[idx] + "\n")
                fileout.write(dic[dickeys[idx]][int(item[0]) - 1: int(item[1])] + "\n")
            fileout.close()
        if len(unreliable_regions) == 1 and lens > unreliable_regions[0][1]:
            fileout = open(dir_output + str(unreliable_regions[0][1]) + "-" + str(lens - 1) + ".reliable", 'w')
            for idx in range(len(dickeys)):
                fileout.write(dickeys[idx] + "\n")
                fileout.write(dic[dickeys[idx]][int(unreliable_regions[0][1]): int(lens)] + "\n")
            fileout.close()
        elif len(unreliable_regions) > 1:
            for i in range(len(unreliable_regions) - 1):
                fileout = open(dir_output + str(unreliable_regions[i][1]) + "-" + str(int(unreliable_regions[i + 1][0]) - 2) + ".reliable", 'w')
                for idx in range(len(dickeys)):
                    fileout.write(dickeys[idx] + "\n")
                    fileout.write(dic[dickeys[idx]][int(unreliable_regions[i][1]): int(unreliable_regions[i + 1][0]) - 1] + "\n")
                fileout.close()
            if unreliable_regions[len(unreliable_regions) -1][1] < lens:
                fileout = open(dir_output + str(unreliable_regions[len(unreliable_regions) -1][1]) + "-" + str(lens - 1) + ".reliable", 'w')
                for idx in range(len(dickeys)):
                    fileout.write(dickeys[idx] + "\n")
                    fileout.write(dic[dickeys[idx]][int(unreliable_regions[len(unreliable_regions) -1][1]): int(lens)] + "\n")
                fileout.close()
