#!/usr/bin/env python

import os
import subprocess
import sys
import math
from joblib import load
sys.path.append("./utils/")
import time
from pandas.core.frame import DataFrame
from utils import Refresh
from utils import Move
from prepare_features_4_classifier_1 import getFeatures4Classifier1
from classifier_c_p_np_aln import AlteredPnp
from calculate_column_scores import calculateColScore
from classifier_realign_strategy import getRealignStrategy
from classifier_region_min_length import getRegionsLength
from seperate_regions import seperateCategory1Regions
from seperate_regions import seperateCategory2Regions
from do_realign import doRealignDir
from do_realign import combineFiles

quickprobs =  "./realign/QuickProbs/bin/quickprobs "
sigma = 1.2
beta = 0.0
threshold = 2.0

realign_normal = "./realign/QuickProbs/bin/quickprobs "
realign_short = "./realign/QuickProbs/bin/quickprobs "

dir_output = "./tmp/seperate_regions/"
col_score = []
output_file = "result.msa"
killed_stage = 0

if __name__ == "__main__":
    start_time = time.time()
    Refresh()
    seq_file = sys.argv[1]
    if len(sys.argv) >= 2:
            output_file = sys.argv[2]
    test_list, prepare_data_1, avg_PID, sd_PID, factor =  getFeatures4Classifier1(seq_file)
    print("[ELAPSED TIME] Preparing data for \"Classifier 1\" takes %.3f sec."%(prepare_data_1 - start_time))
    
    result_real_output, base_msa_time, class1_time, killed_stage = AlteredPnp(test_list, killed_stage, prepare_data_1, seq_file)
    print("[ELAPSED TIME] Get base MSA spends %.3f sec."%(base_msa_time - class1_time))
    
    prepare_data_2, col_score, un_sp, len_seqs, len_family, sd_un_sp, peak_length_ratio = calculateColScore(result_real_output)
    print("[ELAPSED TIME] Preparing data for \"Classifier 3\" takes %.3f sec."%(prepare_data_2 - base_msa_time)) 
    class_region = getRealignStrategy(peak_length_ratio, avg_PID, sd_un_sp, un_sp)
    if int(class_region) == 0:
        print("[MAIN STEP] Choose to run \"Realign Credible Regions(RCR)\" module!")
    else:
        print("[MAIN STEP] Choose to run \"Realign Incredible Regions(RIR)\" module!")

    class3_time = time.time()
    print("[ELAPSED TIME] \"Classifier 3\" takes %.3f sec." % (class3_time - prepare_data_2))
    
    class2_time = 0
    if int(class_region) == 1:
        class_lens = getRegionsLength(len_seqs, len_family, avg_PID, sd_PID, un_sp)
        classifier2_time = time.time()
        print("[ELAPSED TIME] \"Classifier 2\" takes %.3f sec."%(classifier2_time - class3_time))
        killed_stage = seperateCategory1Regions(seq_file, col_score, sigma, beta, class_lens, result_real_output, dir_output, output_file, killed_stage)
        class2_time = time.time()
        print("[ELAPSED TIME] RIR spends %.3f sec."%(class2_time - class3_time))
    else:
        killed_stage = seperateCategory2Regions(seq_file, col_score, threshold, result_real_output, dir_output, output_file, killed_stage)
        class2_time = time.time()
        print("[ELAPSED TIME] RCR spends %.3f sec."%(class2_time - class3_time))
    
    if killed_stage != 4:
        print("[MAIN STEP] Realign !!!")
        doRealignDir(seq_file, dir_output, realign_normal, realign_short, class_region, factor)
        realign_time = time.time()
        print("[ELAPSED TIME] Realigments spend %.3f sec."%(realign_time - class2_time))
        print("[MAIN STEP] Combination !!!")
        
        combineFiles(seq_file, dir_output, output_file)
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
