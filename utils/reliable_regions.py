#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:19:32 2019

@author: mmkuang
"""
import os

quick = "./realign/quickprobs "

def Quickprobs(seq_file, dir_output):
    os.system(quick + seq_file + " > " + dir_output + "quickprobs.txt")
    filein = open(dir_output + "quickprobs.txt", 'r')
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
    with open(dir_output + "quickprobs.txt", 'w') as fileout:
        for dic_ in dickeys:
            fileout.write(dic_ + "\n")
            fileout.write(dic[dic_] + "\n")
    os.system("mv " + dir_output + "quickprobs.txt " + dir_output + "0-" + str(lens) + ".reliable")

def GetReliableRegions(col_score, threshold, class_lens_, seq_file):
    divide_lens = 10
    last_col = len(col_score) - 1
    reliable_regions = []
    tmp_1 = 0
    tmp_2 = 0
    tmp_head = 0
    for item in range(len(col_score)):
        if float(col_score[item]) >= threshold and tmp_1 == 0:
            tmp_head = int(item) + 1
            tmp_1 = 1
        elif float(col_score[item]) >= threshold and tmp_1 == 1 and tmp_2 == 0:
            tmp_2 = 1
        elif float(col_score[item]) >= threshold and tmp_1 == 1 and tmp_2 == 1:
            if item == last_col:
                if int(item) - int(tmp_head) > divide_lens:
                    reliable_regions.append([int(tmp_head), int(item)])
            else:
                continue
        elif float(col_score[item]) < threshold and tmp_1 == 1 and tmp_2 == 1:
            if int(item) - int(tmp_head) > divide_lens:
                reliable_regions.append([int(tmp_head), int(item)])
            tmp_1 = 0
            tmp_2 = 0
            tmp_head = 0
        else:
            tmp_1 = 0
            tmp_2 = 0
            tmp_head = 0
    return reliable_regions

def seperateReliableRegions(reliable_regions, real_output, dir_output):
    file_context = real_output.split("\n")
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
    if len(reliable_regions) == 0:
        with open(dir_output + "0-" + str(lens - 1) + ".reliable", 'w') as filein:
            for line in file_context:
                filein.write(line + "\n")
    else:
        if int(reliable_regions[0][0]) > 1:
            fileout = open(dir_output + "0-" + str(int(reliable_regions[0][0]) - 2) + ".reliable", 'w')
            for idx in range(len(dickeys)):
                fileout.write(dickeys[idx] + "\n")
                fileout.write(dic[dickeys[idx]][0: int(reliable_regions[0][0]) - 1] + "\n")
            fileout.close()
        for item in reliable_regions:
            fileout = open(dir_output + str(int(item[0]) - 1) + "-" + str(int(item[1]) - 1) + ".unreliable", 'w')
            for idx in range(len(dickeys)):
                fileout.write(dickeys[idx] + "\n")
                fileout.write(dic[dickeys[idx]][int(item[0]) - 1: int(item[1])] + "\n")
            fileout.close()
        if len(reliable_regions) == 1 and lens > reliable_regions[0][1]:
            fileout = open(dir_output + str(reliable_regions[0][1]) + "-" + str(lens - 1) + ".reliable", 'w')
            for idx in range(len(dickeys)):
                fileout.write(dickeys[idx] + "\n")
                fileout.write(dic[dickeys[idx]][int(reliable_regions[0][1]): int(lens)] + "\n")
            fileout.close()
        elif len(reliable_regions) > 1:
            for i in range(len(reliable_regions) - 1):
                fileout = open(dir_output + str(reliable_regions[i][1]) + "-" + str(int(reliable_regions[i + 1][0]) - 2) + ".reliable", 'w')
                for idx in range(len(dickeys)):
                    fileout.write(dickeys[idx] + "\n")
                    fileout.write(dic[dickeys[idx]][int(reliable_regions[i][1]): int(reliable_regions[i + 1][0]) - 1] + "\n")
                fileout.close()
            if reliable_regions[len(reliable_regions) -1][1] < lens:
                fileout = open(dir_output + str(reliable_regions[len(reliable_regions) -1][1]) + "-" + str(lens - 1) + ".reliable", 'w')
                for idx in range(len(dickeys)):
                    fileout.write(dickeys[idx] + "\n")
                    fileout.write(dic[dickeys[idx]][int(reliable_regions[len(reliable_regions) -1][1]): int(lens)] + "\n")
                fileout.close()
