#!/usr/bin/env python

import os
import sys
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
from get_TC_SP import getTC_SP
qscore = "./qscore/qscore "
pnp_getpid_path = "./PnpProbs/alter_pnpprobs -G "
pnp_getmsa_path = "./PnpProbs/alter_pnpprobs -p "
real_pnp_output = " -o ./tmp/head_ret.msa"
real_output = "./tmp/head_ret.msa"
pid_path = "./tmp/tmp_pid.txt"
alternative_msa_path = " ./tmp/alternative_msa/ "
tmp_tail_path = "./tmp/tail_ret.seq"
calc_col_score_prefix = " ./tmp/calc_col_score/ -d "
dir_output = "./tmp/seperate_regions/"
col_score = "./tmp/calc_col_score/_col_col.scr"
sigma = 0.7
beta = -1
theta = 0.4
threshold = 0.8
output_file = "result.msa"
killed_stage = 0
model_ = "./classifier/model/abc.joblib"
para_ = "./classifier/model/para.txt"

def getPID(seq_file):
    global killed_stage
    os.system(pnp_getpid_path + seq_file)
    ret_list = []
    if not os.path.exists(pid_path):
        killed_stage = 1
        return [ret_list]
    else:
        if not os.path.getsize(pid_path):
            killed_stage = 1
            return [ret_list]
    para = [0.0 for i in range(6)]
    with open(para_, 'r') as para_in:
        para_context = para_in.read().splitlines()
        for i in range(len(para)):
            para[i] = float(para_context[i])

    with open(pid_path, 'r') as filein:
        file_context = filein.read().splitlines()[0]
        tmp_list = file_context.split("\t")
        for i in range(len(tmp_list)):
            if para[i * 2] - para[i * 2 + 1] > 0:
                ret_list.append((float(tmp_list[i]) - para[i * 2 + 1]) / (para[i * 2] - para[i * 2 + 1]))
            else:
                ret_list.append(float(tmp_list[i]))
    print("Already get classification data.")
    return [ret_list]

def TestClassifier(test_list):
    global killed_stage
    if killed_stage == 1:
        return 0
    clf = load(model_)
    result = []
    result.append(2)
    result = clf.predict(test_list)
    if int(result[0]) == 2:
        return 0
    if int(result[0]) == 0:
        print("Adapt to Progressive PnpProbs.")
    else:
        print("Adapt to non-Progressive PnpProbs.")
    return int(result[0])

def getMSA(class_, seq_file):
    global killed_stage
    print("MSA process is begining ...")
    os.system(pnp_getmsa_path + str(class_) + " " + seq_file + real_pnp_output)
    if not os.path.exists(real_output):
        killed_stage = 2
        return
    else:
        if not os.path.getsize(real_output):
            killed_stage = 2
            return
    processingHead_MSA(real_output)
    print("MSA process ended.")

def getAlternativeMSA(class_, seq_file):
    global killed_stage
    if killed_stage == 2:
        return
    tmp_file_num = 0
    getTail(seq_file, tmp_tail_path)
    print("Alternative MSA processes is begining ...")
    os.system(pnp_getmsa_path + str(1 - int(class_)) + " " + seq_file + " -o ./tmp/alternative_msa/1.msa")
    os.system(pnp_getmsa_path + str(class_) + " " + tmp_tail_path + " -o ./tmp/alternative_msa/2.msa")
    os.system(pnp_getmsa_path + str(1 - int(class_)) + " " + tmp_tail_path + " -o ./tmp/alternative_msa/3.msa")
    if not os.path.exists("./tmp/alternative_msa/1.msa"):
        tmp_file_num += 1
    else:
        if not os.path.getsize("./tmp/alternative_msa/1.msa"):
            tmp_file_num += 1
            os.system("rm ./tmp/alternative_msa/1.msa")
        else:
            processingHead_MSA("./tmp/alternative_msa/1.msa")

    if not os.path.exists("./tmp/alternative_msa/2.msa"):
        tmp_file_num += 1
    else:
        if not os.path.getsize("./tmp/alternative_msa/2.msa"):
            tmp_file_num += 1
            os.system("rm ./tmp/alternative_msa/2.msa")
        else:
            reverseTail_MSA("./tmp/alternative_msa/2.msa")

    if not os.path.exists("./tmp/alternative_msa/3.msa"):
        tmp_file_num += 1
    else:
        if not os.path.getsize("./tmp/alternative_msa/3.msa"):
            tmp_file_num += 1
            os.system("rm ./tmp/alternative_msa/3.msa")
        else:
            reverseTail_MSA("./tmp/alternative_msa/3.msa")

    print("Alternative MSA processes ended.")
    if tmp_file_num < 3:
        print("Calculating HoT Column Scores...")
        os.system("./calc_col_score/calc_col_score " + real_output + calc_col_score_prefix + alternative_msa_path)
        print("Calculated HoT Column Scores.")
    else:
        killed_stage = 3


def Refresh():
    os.system("rm -rf ./tmp")
    os.system("mkdir ./tmp")
    os.system("mkdir ./tmp/calc_col_score")
    os.system("mkdir ./tmp/alternative_msa")
    os.system("mkdir ./tmp/seperate_regions")

def seperateRegions(seq_file, col_score, sigma, beta):
    global killed_stage
    # if not (os.path.exists("./tmp/alternative_msa/1.msa") and os.path.exists("./tmp/alternative_msa/2.msa") and os.path.exists("./tmp/alternative_msa/3.msa")):
    #     killed_stage = 3
    if killed_stage != 2:
        if killed_stage != 3:
            if not os.path.exists(col_score):
                killed_stage = 4
                return
            else:
                if not os.path.getsize(col_score):
                    killed_stage = 4
                    return
            print("Seperating Regions...")
            unreliable_regions = getUnreliableRegions(sigma, beta, theta, threshold, col_score, seq_file, real_output)
            seperateUnreliableRegions(unreliable_regions, real_output, dir_output)
            print("Seperated Regions.")
        else:
            killed_stage = 4
    else:
        killed_stage = 4
        Align_ClustalW2(seq_file)

def Align_ClustalW2(seq_file):
    print("Using ClustalW2 to continue MSA process ...")
    if os.path.exists("./tmp/clustalw"):
        os.system("rm -rf ./tmp/clustalw")
    clustalw2 = "./clustalw/clustalw2 "
    os.system("mkdir ./tmp/clustalw")
    os.system("cp " + seq_file + " ./tmp/clustalw/tmp.seq")
    os.system(clustalw2 + " ./tmp/clustalw/tmp.seq > ./tmp/clustalw/tmp.log")
    tmp_file_in = "./tmp/clustalw/tmp.aln"

    if not os.path.exists(tmp_file_in):
        Move(seq_file, output_file)
    else:
        if not os.path.getsize(tmp_file_in):
            Move(seq_file, output_file)
        else:
            tmp_file_out = "./tmp/clustalw/tmp.msa"
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
                if tmp_list[0] not in dic.keys():
                    dic[tmp_list[0]] = tmp_list[1]
                else:
                    dic[tmp_list[0]] += tmp_list[1]

    dickeys = sorted(dic.keys())
    fileout = open(file_out, 'w')
    for idx in range(len(dickeys)):
        fileout.write(">" + dickeys[idx] + "\n")
        fileout.write(dic[dickeys[idx]] + "\n")
    fileout.close()

def CheckFile(output_file, class_):
    need_remake = False
    dic = {}
    has_key = False
    value = ""
    key = ""
    filein = open(output_file, 'r')
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
    dickeys = sorted(dic.keys())
    for idx in range(len(dickeys)):
        if len(dic[dickeys[idx]]) != len(dic[dickeys[0]]):
            need_remake = True
    if need_remake == True:
        os.system(pnp_getmsa_path + str(class_) + " " + seq_file + " -o " + output_file)

def Move(file_in, file_out):
    os.system("cp "+ file_in + " " + file_out)

if __name__ == "__main__":
    Refresh()
    seq_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    test_list = getPID(seq_file)
    class_ = TestClassifier(test_list)
    getMSA(class_, seq_file)
    # getAlternativeMSA(class_, seq_file)
    detect_unreliable_regions(real_output, col_score)
    seperateRegions(seq_file, col_score, sigma, beta)
    if killed_stage != 4:
        print("Realign !!!")
        do_Realign_Dir(dir_output, class_, pnp_getmsa_path)
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
                Align_ClustalW2(seq_file)
            else:
                if not os.path.getsize(output_file):
                    Align_ClustalW2(seq_file)
    killed_stage = 0
    #Move(real_output, output_file)
    if not os.path.getsize(output_file):
        print("Result is Empty ?")
        Align_ClustalW2(seq_file)
