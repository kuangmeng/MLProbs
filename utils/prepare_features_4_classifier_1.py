#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:19:32 2019

@author: mmkuang
"""
import os
import subprocess
import sys
import time
pnp_getpid_path = "./baseMSA/PnpProbs/alter_pnpprobs -G "

para_ = "./classifier/model/branch/para.txt"

def readPID(pid_out):
    file_context = pid_out.split("\t")
    if len(file_context) >= 7:
        return file_context[1], file_context[4], file_context[5], [file_context[0], file_context[2], file_context[3]], file_context[6]
    else:
        return 0, 0, 0,  [0, 0, 0], 0

def getFeatures4Classifier1(seq_file):
    rc, pid_out = subprocess.getstatusoutput(pnp_getpid_path + seq_file)
    prepare_data_1 = time.time()
    sd_PID, avg_sp, peak_length_ratio, tmp_list1, factor = readPID(pid_out)
    avg_PID = float(tmp_list1[0])
    tmp_list = []
    for i in tmp_list1:
        tmp_list.append(i)
    tmp_list.append(avg_sp)
    tmp_list.append(peak_length_ratio)
    ret_list = []
    para = []
    with open(para_, 'r') as para_in:
        para_context = para_in.read().splitlines()
        for i in range(len(para_context)):
            para.append(float(para_context[i]))
    for i in range(len(tmp_list)):
        ret_list.append((float(tmp_list[i]) - para[i * 2 + 1])/ (para[i * 2] - para[i * 2 + 1]))
    print("[MAIN STEP] Already get classification data.")
    return [ret_list], prepare_data_1, avg_PID, sd_PID, factor
