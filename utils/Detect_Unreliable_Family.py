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

    if need_realign == True:
        Align_ClustalW2(seq_file, output_file)

    return need_realign


def Align_ClustalW2(seq_file, output_file):
    print("Low Col Score ! Using ClustalW2 to Realign ...")
    if os.path.exists("./tmp/re_clustalw"):
        os.system("rm -rf ./tmp/re_clustalw")
    clustalw2 = "./clustalw/clustalw2 "
    os.system("mkdir ./tmp/re_clustalw")
    os.system("cp " + seq_file + " ./tmp/re_clustalw/tmp.seq")
    os.system(clustalw2 + " ./tmp/re_clustalw/tmp.seq > ./tmp/re_clustalw/tmp.log")
    tmp_file_in = "./tmp/re_clustalw/tmp.aln"

    if  os.path.exists(tmp_file_in) and os.path.getsize(tmp_file_in):
        tmp_file_out = "./tmp/re_clustalw/tmp.msa"
        Calc_FASTA(tmp_file_in, tmp_file_out)
        Move(tmp_file_out, output_file)
    print("ClustalW2 finished.")

def Calc_FASTA(file_in, file_out):
    dic = {}
    with open(file_in, 'r') as filein:
        file_context = filein.read().splitlines()
        for item in file_context:
            tmp_list = item.split()
            if len(tmp_list) == 2:
                if ("*" in tmp_list[0]) or (":" in tmp_list[0]):
                    continue
                elif tmp_list[0] not in dic.keys():
                    dic[tmp_list[0]] = tmp_list[1]
                else:
                    dic[tmp_list[0]] += tmp_list[1]

    dickeys = sorted(dic.keys())
    fileout = open(file_out, 'w')
    for idx in range(len(dickeys)):
        fileout.write(">" + dickeys[idx] + "\n")
        fileout.write(dic[dickeys[idx]] + "\n")
    fileout.close()


def Move(file_in, file_out):
    os.system("cp "+ file_in + " " + file_out)
