#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:19:32 2019

@author: mmkuang
"""
from joblib import load
import time
import subprocess

model_ = "./classifier/model/branch/randomforest.joblib"
para_ = "./classifier/model/branch/para.txt"
pnp_getmsa_path = "./baseMSA/C_P_NP_Aln/c_p_np_aln -p "
quickprobs =  "./realign/QuickProbs/bin/quickprobs "

def testClassifier(test_list, killed_stage):
    global model_
    global para_
    if killed_stage == 1:
        return 0
    clf = load(model_)
    result = clf.predict(test_list)
    if int(result[0]) >= 2 or int(result[0]) < 0:
        return 0
    if int(result[0]) == 0:
        print("[MAIN STEP] Adapt to Progressive Strategy.")
    else:
        print("[MAIN STEP] Adapt to non-Progressive Strategy.")
    return int(result[0])

def getMSA(class_, seq_file, killed_stage):
    print("[MAIN STEP] MSA process is begining ...")
    status = 0
    result_real_output = ""
    if class_ < 2:
        status, result_real_output = subprocess.getstatusoutput(pnp_getmsa_path + str(class_) + " " + seq_file)
    else:
        status, result_real_output = subprocess.getstatusoutput(quickprobs + " " + seq_file)
    if status != 0:
        killed_stage  = 2
    print("[MAIN STEP] MSA process ended.")
    return result_real_output, killed_stage

def AlteredPnp(test_list, killed_stage, prepare_data_1, seq_file):
    class_ = testClassifier(test_list, killed_stage)
    class1_time = time.time()
    print("[ELAPSED TIME] \"Classifier 1\" takes %.3f sec."%(class1_time - prepare_data_1))
    result_real_output, killed_stage = getMSA(class_, seq_file, killed_stage)
    base_msa_time = time.time() 
    return result_real_output, base_msa_time, class1_time, killed_stage
