#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 17:21:35 2019

@author: mmkuang
"""

def DeleteN(seq_file):
    filein = open(seq_file, 'r')
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

    return dic

def getTail(seq_file, save_path_tail):
    dic = DeleteN(seq_file)
    dickeys = sorted(dic.keys())
    fileout = open(save_path_tail, 'w')
    for idx in range(len(dickeys)):
        fileout.write(dickeys[idx] + "\n")
        fileout.write(dic[dickeys[idx]][::-1] + "\n")
    fileout.close()
    
        