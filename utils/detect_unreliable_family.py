#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 15:42:00 2019

@author: mmkuang
"""
import os

def Detect_Unreliable(theta, threshold, col_score, seq_file, output_file):
    num_score = []
    need_realign = False
    with open(col_score, 'r') as filein:
        file_context = filein.read().splitlines()
        for i in range(1, len(file_context) - 1):
            num_score.append(file_context[i].split())
    lens = len(num_score)
    tmp_unreliable = 0
    for item in num_score:
        if float(item[1]) <= theta:
            tmp_unreliable += 1
    if float(tmp_unreliable) / float(lens) >= threshold:
        need_realign = True
        
    return need_realign


def Move(file_in, file_out):
    os.system("cp "+ file_in + " " + file_out)
