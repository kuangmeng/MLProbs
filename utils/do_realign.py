#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 18:49:18 2019

@author: mmkuang
"""
import os
import re
from Detect_Unreliable_Family import Align_ClustalW2

def Realign_file_list(dir_output):
    file_list = os.listdir(dir_output)
    realign_file_list = []
    for file in file_list:
        if os.path.splitext(file)[-1][1:] == "unreliable" and file[0] != '.':
            realign_file_list.append(dir_output + file)
    return realign_file_list

def perprocess(real_output, real_out, tmp_array):
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
    dickeys = sorted(dic.keys())
    with open(real_out, 'w') as fileout:
        for item in dickeys:
            if bool(re.search('[A-Z]', dic[item])):
                fileout.write(item + "\n")
                fileout.write(dic[item].replace("-", "").replace(".", "") + "\n")
            else:
                tmp_array.append(item)

def do_Realign(class_, quick, realign_short, which_part, realign_file):
    ret_name = ""
    tmp_path = realign_file.split("/")
    for i in range(len(tmp_path) - 1):
        ret_name += tmp_path[i] + "/"
    ret_name += os.path.splitext(tmp_path[len(tmp_path) - 1])[0] + ".reliable"
    if not os.path.exists("./tmp/qp_tmp"):
        os.system("mkdir ./tmp/qp_tmp/")
    tmp_file = "./tmp/qp_tmp/" + os.path.splitext(tmp_path[len(tmp_path) - 1])[0] + ".unreliable"
    tmp_array = []
    perprocess(realign_file, tmp_file, tmp_array)
    #os.system(pnp_path + str(class_) + " " + realign_file + " -o " + ret_name)
    if int(which_part) == 0:
        os.system(realign_short + " " + tmp_file + " > " + ret_name)
    else:
        os.system(quick + " " + tmp_file + " > " + ret_name)
    if not os.path.exists(ret_name):
        os.system("cp " +  realign_file + "  " + ret_name)
    else:
        if not os.path.getsize(ret_name):
            os.system("cp " + realign_file  + "  " + ret_name)

    add_PerProcess(ret_name, tmp_array)

def add_PerProcess(ret_name, tmp_save):
    filein = open(ret_name, 'r')
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
    dickeys = sorted(dic.keys())
    lens = len(dic[dickeys[0]])
    with open(ret_name, 'w') as fileout:
        for item in dickeys:
            fileout.write(item + "\n")
            fileout.write(dic[item] + "\n")
        for itm in tmp_save:
            fileout.write(itm + "\n")
            fileout.write("-"*lens + "\n")

def do_Realign_Dir(dir_output, class_, quick, realign_short, which_part):
    realign_file_list = Realign_file_list(dir_output)
    for file in realign_file_list:
        do_Realign(class_, quick, realign_short, which_part, file)

def getFile_Len(filename):
    seq_file_lens = 0
    with open(filename, 'r') as seq_file_in:
        seq_file_context = seq_file_in.read().splitlines()
        for item in seq_file_context:
            if item.strip()[0:1] == ">":
                seq_file_lens += 1
    return seq_file_lens

def Combination_Files(seq_file, dir_output, output_file):
    seq_file_lens = getFile_Len(seq_file)
    file_list = os.listdir(dir_output)
    need_combination = []
    for file in file_list:
        if os.path.splitext(file)[-1][1:] == "reliable" and file[0] != '.':
            need_combination.append(dir_output + file)
    tmp_num_arr = []
    for itt in need_combination:
        tmp_num_arr.append(int(itt.split("/")[-1].split("-")[0]))
    tmp_num_arr = sorted(tmp_num_arr)
    need_combination_files = []
    for num_ in tmp_num_arr:
        for n_file in need_combination:
            if str(num_) == n_file.split("/")[-1].split("-")[0]:
                need_combination_files.append(n_file)
    if len(need_combination) != len(need_combination_files):
        print("ERROR: file length")
        return

    dic = {}
    has_key = False
    value = ""
    key = ""
    tmp_file_name = need_combination_files[0]
    tmp_file_lens = getFile_Len(tmp_file_name)
    if (not os.path.getsize(tmp_file_name)) or tmp_file_lens != seq_file_lens:
        tmp_path_list = tmp_file_name.split(".")
        tmp_file_name = "." + tmp_path_list[1] + ".unreliable"
        # Align_ClustalW2(tmp_file_name_, tmp_file_name)
        # if (not os.path.getsize(tmp_file_name)) or tmp_file_lens != seq_file_lens:
        #     tmp_path_list = tmp_file_name.split(".")
        #     tmp_file_name = "." + tmp_path_list[1] + ".unreliable"
        print("Fixed: No sequences read Error !")
    filein = open(tmp_file_name, 'r')
    file_context = filein.read().splitlines()
    filein.close()
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
    if len(need_combination_files) > 1:
        for file_i in range(1, len(need_combination_files)):
            file_name = need_combination_files[file_i]
            file_lens = getFile_Len(file_name)
            if (not os.path.getsize(file_name)) or file_lens != seq_file_lens:
                path_list = file_name.split(".")
                file_name = "." + path_list[1] + ".unreliable"
                # Align_ClustalW2(file_name_, file_name)
                # if (not os.path.getsize(file_name)) or tmp_file_lens != seq_file_lens:
                #     path_list = file_name.split(".")
                #     file_name = "." + path_list[1] + ".unreliable"
                print("Fixed: No sequences read Error !")
            file_in = open(file_name)
            filein_context = file_in.read().splitlines()
            file_in.close()
            tmp_key = ""
            tmp_value = ""
            tmp_haskey = False
            for t in range(len(filein_context)):
                if filein_context[t][0:1] == ">":
                    if tmp_haskey == True:
                        dic[tmp_key] += tmp_value
                        tmp_value = ""
                        tmp_key = ""
                        tmp_haskey = False
                    tmp_haskey = True
                    tmp_key = filein_context[t]
                elif tmp_haskey == True:
                    tmp_value = tmp_value.replace("\r","") + filein_context[t].replace("\r","")
            dic[tmp_key] += tmp_value

    dickeys = sorted(dic.keys())
    fileout = open(output_file, 'w')

    for idx in range(len(dickeys)):
        fileout.write(dickeys[idx] + "\n")
        fileout.write(dic[dickeys[idx]] + "\n")
    fileout.close()
