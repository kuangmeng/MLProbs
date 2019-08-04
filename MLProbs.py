#!/usr/bin/env python

import os
import sys
import math
from joblib import load
sys.path.append("./utils/")
sys.path.append("./classifier/")
from do_realign import do_Realign_Dir
from do_realign import Combination_Files
from unreliable_regions import getUnreliableRegions
from unreliable_regions import seperateUnreliableRegions
from preprocessing_seq_file import getTail
from postprocessing_msa_file import processingHead_MSA
from postprocessing_msa_file import reverseTail_MSA
from Detect_Unreliable_Regions import detect_unreliable_regions
from reliable_regions import Quickprobs
from reliable_regions import GetReliableRegions
from reliable_regions import seperateReliableRegions
quickprobs = "./realign/quickprobs "
pnp_getpid_path = "./baseMSA/pnpprobs/alter_pnpprobs -G "
pnp_getmsa_path = "./baseMSA/pnpprobs/alter_pnpprobs -p "
real_pnp_output = " -o ./tmp/head_ret.msa"
real_output = "./tmp/head_ret.msa"
pid_path = "./tmp/tmp_pid.txt"
tmp_tail_path = "./tmp/tail_ret.seq"
calc_col_score_prefix = " ./tmp/calc_col_score/ -d "
dir_output = "./tmp/seperate_regions/"
col_score = "./tmp/calc_col_score/_col_col.scr"
sigma = 1.0
beta = 0.0
theta = 0.4
threshold = 0
output_file = "result.msa"
killed_stage = 0
model_ = "./classifier/model/branch/abc.joblib"
para_ = "./classifier/model/branch/para.txt"
which_part = 0
model_lens = "./classifier/model/seq_lens/randomforest.joblib"
para_lens = "./classifier/model/seq_lens/para.txt"
model_lens_short = "./classifier/model/seq_lens_short/randomforest.joblib"
para_lens_short = "./classifier/model/seq_lens_short/para.txt"

model_region = "./classifier/model/regions/randomforest.joblib"
para_region = "./classifier/model/regions/para.txt"
model_region_short = "./classifier/model/regions_short/randomforest.joblib"
para_region_short = "./classifier/model/regions_short/para.txt"

mafft = "mafft"

class_region = 0

avg_PID = 0.0
len_seqs = 0
un_sp = 0.0
sd_PID = 0.0
len_family = 0
sd_un_sp = 0.0
peak_length_ratio = 0.0


def ReadFile():
    with open("./tmp/tmp_pid.txt", 'r') as filein:
        file_context = filein.read().splitlines()[0].split("\t")
    if len(file_context) >= 4:
        return file_context[1], [file_context[0], file_context[2], file_context[3]]
    else:
        return 0, [0, 0, 0]

def Calculate(len_family, len_seq):
    if len_family == 0 or len_seq == 0:
        return 0, 0, 0
    dic = {}
    dir_ = "./tmp/tmp_pid/"
    file_list = os.listdir(dir_)
    for file in file_list:
        if file[0] != '.':
            with open(dir_ + file) as filein:
                tmp_list = filein.read().split("\t")[:-1]
                for idx in range(0, len(tmp_list), 2):
                    if tmp_list[idx] in dic.keys():
                        dic[tmp_list[idx]] += float(tmp_list[idx + 1])
                    else:
                        dic[tmp_list[idx]] = float(tmp_list[idx + 1])
    divided = ((len_family - 1) * len_family) / 2
    avg_sp = 0.0
    dickeys = sorted(dic.keys())
    for i in range(len(dickeys)):
        if int(dickeys[i]) < len_seq:
            dic[dickeys[i]] /= divided
            avg_sp += dic[dickeys[i]]
    avg_sp /= len_seq
    sd_sp = 0.0
    for i in range(len(dickeys)):
        if int(dickeys[i]) < len_seq:
            sd_sp += ( dic[dickeys[i]] - avg_sp ) ** 2
    sd_sp /= len_seq
    sd_sp = math.sqrt(sd_sp)
    peak_length_ratio = 0.0
    for i in range(len(dickeys)):
        if int(dickeys[i]) < len_seq:
            if dic[dickeys[i]] >= 1:
                peak_length_ratio += 1
    peak_length_ratio /= len_seq
    if sd_sp > 10:
        sd_sp = 0
    if avg_sp > 10:
        avg_sp = 0
    if peak_length_ratio > 10:
        peak_length_ratio = 0
    return [avg_sp, peak_length_ratio]

def getPID(seq_file, which_part):
    global avg_PID
    global sd_PID
    global killed_stage
    os.system(pnp_getpid_path + seq_file)
    sd_PID, tmp_list1 = ReadFile()
    avg_PID = float(tmp_list1[0])
    if which_part == 0:
        return [avg_PID]
    tmp_list2 = Calculate(int(tmp_list1[1]), int(tmp_list1[2]))
    tmp_list = []
    for i in tmp_list1:
        tmp_list.append(i)
    for i in tmp_list2:
        tmp_list.append(i)
    ret_list = []
    if not os.path.exists(pid_path):
        killed_stage = 1
        return [ret_list]
    else:
        if not os.path.getsize(pid_path):
            killed_stage = 1
            return [ret_list]
    para = []
    with open(para_, 'r') as para_in:
        para_context = para_in.read().splitlines()
        for i in range(len(para_context)):
            para.append(float(para_context[i]))
    for i in range(len(tmp_list)):
        ret_list.append((float(tmp_list[i]) - para[i * 2 + 1])/ (para[i * 2] - para[i * 2 + 1]))
    print("Already get classification data.")
    return [ret_list]

def TestClassifier(test_list):
    global killed_stage
    if killed_stage == 1:
        return 0
    clf = load(model_)
    result = clf.predict(test_list)
    if int(result[0]) >= 2 or int(result[0]) < 0:
        return 0
    if int(result[0]) == 0:
        print("Adapt to Progressive PnpProbs.")
    else:
        print("Adapt to non-Progressive PnpProbs.")
    return int(result[0])

def getMSA(class_, seq_file):
    global killed_stage
    print("MSA process is begining ...")
    if class_ < 2:
        os.system(pnp_getmsa_path + str(class_) + " " + seq_file + real_pnp_output)
    else:
        os.system(mafft + " " + seq_file + " > " + real_output)
    if not os.path.exists(real_output):
        killed_stage = 2
        return
    else:
        if not os.path.getsize(real_output):
            killed_stage = 2
            return
    processingHead_MSA(real_output)
    print("MSA process ended.")

def Refresh():
    os.system("rm -rf ./tmp")
    os.system("mkdir ./tmp")
    os.system("mkdir ./tmp/calc_col_score")
    os.system("mkdir ./tmp/seperate_regions")
    os.system("mkdir ./tmp/tmp_pid")

def seperateRegions(seq_file, col_score, sigma, beta, class_lens):
    global killed_stage
    if killed_stage != 2:
        if killed_stage != 3:
            if not os.path.exists(col_score):
                killed_stage = 4
                return
            else:
                if not os.path.getsize(col_score):
                    killed_stage = 4
                    return
            print("Seperating Unreliable Regions...")
            unreliable_regions = getUnreliableRegions(sigma, beta, theta, threshold, col_score, seq_file, real_output, class_lens)
            seperateUnreliableRegions(unreliable_regions, real_output, dir_output)
            print("Seperated Unreliable Regions.")
        else:
            killed_stage = 4
    else:
        killed_stage = 4
        os.system(quickprobs + " " + seq_file + " > " + output_file)

def SeperateReliableRegions(seq_file, col_score, which_part):
    global killed_stage
    if killed_stage != 2:
        if killed_stage != 3:
            if not os.path.exists(col_score):
                killed_stage = 4
                return
            else:
                if not os.path.getsize(col_score):
                    killed_stage = 4
                    return
            print("Seperating Reliable Regions...")
            if int(which_part) == 1:
                Quickprobs(seq_file, dir_output)
            else:
                reliable_regions = GetReliableRegions(col_score, threshold, 0, seq_file)
                seperateReliableRegions(reliable_regions, real_output, dir_output)
            print("Seperated Reliable Regions.")
        else:
            killed_stage = 4
    else:
        killed_stage = 4
        os.system(quickprobs + " " + seq_file + " > " + output_file)

def Move(file_in, file_out):
    os.system("cp "+ file_in + " " + file_out)

def getRegionsLength(len_seqs, len_family, avg_PID, sd_PID, un_sp, which_part):
    class_lens = 3
    test = [len_seqs, len_family, avg_PID, sd_PID, un_sp]
    para = []
    if int(which_part) == 0:
        para_lens = para_lens_short
        model_lens = model_lens_short
    with open(para_lens, 'r') as filein:
        file_context = filein.read().splitlines()
        for i in range(2 * len(test)):
            para.append(float(file_context[i]))
    real_test = []
    for i in range(len(test)):
        real_test.append((float(test[i]) - para[i * 2 + 1])/ (para[i * 2] - para[i * 2 + 1]))
    clf = load(model_lens)
    class_lens = clf.predict([real_test])[0]
    if class_lens > 3 or class_lens < 0:
        class_lens = 3
    return class_lens

def getRegions_to_Realign(peak_length_ratio, avg_PID, sd_un_sp, un_sp, which_part):
    class_region = 1
    test = [peak_length_ratio, avg_PID, sd_un_sp, un_sp]
    para = []
    if which_part == 0:
        para_region = para_region_short
        model_region = model_region_short
    with open(para_region, 'r') as filein:
        file_context = filein.read().splitlines()
        for i in range(2 * len(test)):
            para.append(float(file_context[i]))
    real_test = []
    for i in range(len(test)):
        real_test.append((float(test[i]) - para[i * 2 + 1])/ (para[i * 2] - para[i * 2 + 1]))
    clf = load(model_region)
    class_region = clf.predict([real_test])[0]
    if class_region > 1 or class_region < 0:
        class_region = 1
    return class_region

if __name__ == "__main__":
    Refresh()
    seq_file = sys.argv[1]
    if len(sys.argv) >= 2:
        if sys.argv[2] == "-s":
            which_part = 0
        elif sys.argv[2] == "-b":
            which_part = 1
        else:
            print("The 2nd option is wrong! Use the default option.")
            which_part = 1
        if len(sys.argv) > 3:
            output_file = sys.argv[3]
    test_list = getPID(seq_file, which_part)
    if which_part == 1:
        class_ = TestClassifier(test_list)
    else:
        class_ = 2
    getMSA(class_, seq_file)
    un_sp, len_seqs, len_family, sd_un_sp, peak_length_ratio = detect_unreliable_regions(real_output, col_score)
    class_region = getRegions_to_Realign(peak_length_ratio, avg_PID, sd_un_sp, un_sp, which_part)
    if int(class_region) == 0:
        print("Choose to Realign Reliable Regions!")
    else:
        print("Choose to Realign Unreliable Regions!")
    if int(class_region) == 1:
        class_lens = getRegionsLength(len_seqs, len_family, avg_PID, sd_PID, un_sp, which_part)
        seperateRegions(seq_file, col_score, sigma, beta, class_lens)
    else:
        SeperateReliableRegions(seq_file, col_score, which_part)

    if killed_stage != 4:
        print("Realign !!!")
        do_Realign_Dir(dir_output, class_, quickprobs)
        print("Combination !!!")
        Combination_Files(seq_file, dir_output, output_file)
        #CheckFile(output_file, class_)
        print("Got the final MSA!")
    else:
        if os.path.exists(real_output) and os.path.getsize(real_output):
            print("Take the original output as Result!")
            Move(real_output, output_file)
        else:
            if not os.path.exists(output_file):
                os.system(quickprobs + " " + seq_file + " > " + output_file)
            else:
                if not os.path.getsize(output_file):
                    os.system(quickprobs + " " + seq_file + " > " + output_file)
    killed_stage = 0
    # Move(real_output, output_file)
    if not os.path.getsize(output_file):
        print("Result is Empty ?")
        os.system(quickprobs + " " + seq_file + " > " + output_file)
