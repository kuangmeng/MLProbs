#!/usr/bin/env python

import os
import subprocess
import sys
import math
from joblib import load
sys.path.append("./utils/")
import time
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

sigma = 1.0
beta = 0.0
theta = 0.4
threshold = 1.0

quickprobs =  "./realign/QuickProbs/bin/quickprobs "
realign_short = "./realign/QuickProbs/bin/quickprobs "

pnp_getpid_path = "./baseMSA/PnpProbs/alter_pnpprobs -G "
pnp_getmsa_path = "./baseMSA/PnpProbs/alter_pnpprobs -p "
dir_output = "./tmp/seperate_regions/"

col_score = []

output_file = "result.msa"

killed_stage = 0

model_ = "./classifier/model/branch/abc.joblib"
para_ = "./classifier/model/branch/para.txt"

which_part = 0

model_lens = "./classifier/model/seq_lens/randomforest.joblib"
para_lens = "./classifier/model/seq_lens/para.txt"

model_region = "./classifier/model/regions/randomforest.joblib"
para_region = "./classifier/model/regions/para.txt"
model_region_short = "./classifier/model/regions_short/randomforest.joblib"
para_region_short = "./classifier/model/regions_short/para.txt"

class_region = 0

avg_PID = 0.0
len_seqs = 0
un_sp = 0.0
sd_PID = 0.0
len_family = 0
sd_un_sp = 0.0
peak_length_ratio = 0.0

def ReadPID(pid_out):
    file_context = pid_out.split("\t")
    if len(file_context) >= 6:
        return file_context[1], file_context[4], file_context[5], [file_context[0], file_context[2], file_context[3]]
    else:
        return 0, 0, 0,  [0, 0, 0]

def getPID(seq_file, which_part):
    global avg_PID
    global sd_PID
    global killed_stage
    rc, pid_out = subprocess.getstatusoutput(pnp_getpid_path + seq_file)
    prepare_data_1 = time.time()
    sd_PID, avg_sp, peak_length_ratio, tmp_list1 = ReadPID(pid_out)
    avg_PID = float(tmp_list1[0])
    if which_part == 0:
        return [avg_PID], prepare_data_1
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
    return [ret_list], prepare_data_1

def TestClassifier(test_list):
    global killed_stage
    if killed_stage == 1:
        return 0
    clf = load(model_)
    result = clf.predict(test_list)
    if int(result[0]) >= 2 or int(result[0]) < 0:
        return 0
    if int(result[0]) == 0:
        print("[MAIN STEP] Adapt to Progressive PnpProbs.")
    else:
        print("[MAIN STEP] Adapt to non-Progressive PnpProbs.")
    return int(result[0])

def getMSA(class_, seq_file):
    global killed_stage
    print("[MAIN STEP] MSA process is begining ...")
    status = 0
    result_real_output = ""
    if class_ < 2:
        status, result_real_output = subprocess.getstatusoutput(pnp_getmsa_path + str(class_) + " " + seq_file)
    else:
        status, result_real_output = subprocess.getstatusoutput(quickprobs + " " + seq_file)
    if status != 0:
        killed_stage  = 2
        return ""
    print("[MAIN STEP] MSA process ended.")
    return result_real_output

def Refresh():
    os.system("rm -rf ./tmp/")
    os.system("mkdir ./tmp/")
    os.system("mkdir ./tmp/seperate_regions")

def seperateRegions(seq_file, col_score, sigma, beta, class_lens, real_output):
    global killed_stage
    if killed_stage != 2:
        if killed_stage != 3:
            print("[MAIN STEP] Seperating Incredible Regions...")
            unreliable_regions = getUnreliableRegions(sigma, beta, theta, threshold, col_score, seq_file, real_output, class_lens)
            seperateUnreliableRegions(unreliable_regions, real_output, dir_output)
            print("[MAIN STEP] Seperated Incredible Regions.")
        else:
            killed_stage = 4
    else:
        killed_stage = 4
        os.system(quickprobs + "  " + seq_file + " > " + output_file)

def SeperateReliableRegions(seq_file, col_score, real_output, which_part):
    global killed_stage
    if killed_stage != 2:
        if killed_stage != 3:
            print("[MAIN STEP] Seperating Credible Regions...")
            if int(which_part) == 0:
                Quickprobs(seq_file, dir_output)
            else:
                reliable_regions = GetReliableRegions(col_score, threshold, 30, seq_file)
                seperateReliableRegions(reliable_regions, real_output, dir_output)
            print("[MAIN STEP] Seperated Credible Regions.")
        else:
            killed_stage = 4
    else:
        killed_stage = 4
        os.system(quickprobs + "  " + seq_file + " > " + output_file)

def Move(file_in, file_out):
    os.system("cp "+ file_in + " " + file_out)

def getRegionsLength(len_seqs, len_family, avg_PID, sd_PID, un_sp, which_part):
    global para_lens
    global model_lens
    class_lens = 3
    test = [len_seqs, len_family, avg_PID, sd_PID, un_sp]
    para = []
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

def getRegions_to_Realign(peak_length_ratio, avg_PID, sd_un_sp, un_sp):
    global para_region
    global model_region
    class_region = 1
    test = [peak_length_ratio, avg_PID, sd_un_sp, un_sp]
    para = []
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


def getRegions_to_Realign_short(peak_length_ratio, sd_un_sp, un_sp):
    global para_region_short
    global model_region_short
    class_region = 1
    test = [peak_length_ratio, sd_un_sp, un_sp]
    para = []
    with open(para_region_short, 'r') as filein:
        file_context = filein.read().splitlines()
        for i in range(2 * len(test)):
            para.append(float(file_context[i]))
    real_test = []
    for i in range(len(test)):
        real_test.append((float(test[i]) - para[i * 2 + 1])/ (para[i * 2] - para[i * 2 + 1]))
    clf = load(model_region_short)
    class_region = clf.predict([real_test])[0]
    if class_region > 1 or class_region < 0:
        class_region = 1
    return class_region


if __name__ == "__main__":
    start_time = time.time()
    Refresh()
    seq_file = sys.argv[1]
    if len(sys.argv) >= 2:
        if sys.argv[2] == "-s":
            which_part = 0
        elif sys.argv[2] == "-b":
            which_part = 1
        else:
            print("[ERROR] The 2nd option is wrong! Use the default option.")
            which_part = 1
        if len(sys.argv) > 3:
            output_file = sys.argv[3]
    class1_time = 0
    if which_part == 1:
        test_list, prepare_data_1 = getPID(seq_file, which_part)
        print("[ELAPSED TIME] Preparing data for \"Classifier 1\" takes %.3f sec."%(prepare_data_1 - start_time))
        class_ = TestClassifier(test_list)
        class1_time = time.time()
        print("[ELAPSED TIME] \"Classifier 1\" takes %.3f sec."%(class1_time - prepare_data_1))
    else:
        class_ = 2
        class1_time = start_time
    result_real_output = getMSA(class_, seq_file)
    base_msa_time = time.time()
    print("[ELAPSED TIME] Get base MSA spends %.3f sec."%(base_msa_time - class1_time))
    prepare_data_2, col_score, un_sp, len_seqs, len_family, sd_un_sp, peak_length_ratio = detect_unreliable_regions(result_real_output)
    print("[ELAPSED TIME] Preparing data for \"Classifier 3\" takes %.3f sec."%(prepare_data_2 - base_msa_time))
    class_region = 0
    if which_part == 1:
        class_region = getRegions_to_Realign(peak_length_ratio, avg_PID, sd_un_sp, un_sp)
    else:
        class_region = getRegions_to_Realign_short(peak_length_ratio, sd_un_sp, un_sp)

    if int(class_region) == 0:
        print("[MAIN STEP] Choose to run \"Realign Credible Regions(RCR)\" module!")
    else:
        print("[MAIN STEP] Choose to run \"Realign Incredible Regions(RIR)\" module!")
    class3_time = time.time()
    print("[ELAPSED TIME] \"Classifier 3\" takes %.3f sec." % (class3_time - prepare_data_2))
    class2_time = 0
    if int(class_region) == 1:
        if int(which_part) == 1:
            class_lens = getRegionsLength(len_seqs, len_family, avg_PID, sd_PID, un_sp, which_part)
            classifier2_time = time.time()
            print("[ELAPSED TIME] \"Classifier 2\" takes %.3f sec."%(classifier2_time - class3_time))
        else:
            class_lens = 30
        seperateRegions(seq_file, col_score, sigma, beta, class_lens, result_real_output)
        class2_time = time.time()
        print("[ELAPSED TIME] RIR spends %.3f sec."%(class2_time - class3_time))
    else:
        SeperateReliableRegions(seq_file, col_score, result_real_output, which_part)
        class2_time = time.time()
        print("[ELAPSED TIME] RCR spends %.3f sec."%(class2_time - class3_time))
    if killed_stage != 4:
        print("[MAIN STEP] Realign !!!")
        do_Realign_Dir(dir_output, class_, quickprobs, realign_short, which_part)
        realign_time = time.time()
        print("[ELAPSED TIME] Realigments spend %.3f sec."%(realign_time - class2_time))
        print("[MAIN STEP] Combination !!!")
        Combination_Files(seq_file, dir_output, output_file)
        print("[MAIN STEP] Got the final MSA!")
        end_time = time.time()
        total_time = end_time - start_time
        print("[ELAPSED TIME] Total Running time: %.3f sec."%(total_time))
    else:
        if not os.path.exists(output_file):
            os.system(quickprobs + " " + seq_file + " > " + output_file)
        else:
            if not os.path.getsize(output_file):
                os.system(quickprobs + " " + seq_file + " > " + output_file)
        end_time = time.time()
        total_time = end_time - start_time
        print("[ELAPSED TIME] Total Running time: %.3f sec."%(total_time))
    killed_stage = 0
    if not os.path.getsize(output_file):
        print("[ERROR] Result is Empty ?")
        os.system(quickprobs + " " + seq_file + " > " + output_file)
        end_time = time.time()
        total_time = end_time - start_time
        print("[ELAPSED TIME] Total Running time: %.3f sec."%(total_time))
