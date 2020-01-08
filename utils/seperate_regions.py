#!/usr/bin/env python
import os
from unreliable_regions import getUnreliableRegions
from unreliable_regions import seperateUnreliableRegions
from reliable_regions import getReliableRegions
from reliable_regions import seperateReliableRegions

quickprobs =  "./realign/QuickProbs/bin/quickprobs "


def seperateCategory1Regions(seq_file, col_score, sigma, beta, class_lens, real_output, dir_output, output_file,  killed_stage):
    global quickprobs
    if killed_stage != 2:
        if killed_stage != 3:
            print("[MAIN STEP] Seperating Incredible Regions...")
            unreliable_regions = getUnreliableRegions(sigma, beta, col_score, seq_file, real_output, class_lens)
            seperateUnreliableRegions(unreliable_regions, real_output, dir_output)
            print("[MAIN STEP] Seperated Incredible Regions.")
        else:
            killed_stage = 4
    else:
        killed_stage = 4
        os.system(quickprobs + "  " + seq_file + " > " + output_file)
    return killed_stage

def seperateCategory2Regions(seq_file, col_score, threshold,  real_output, dir_output, output_file, killed_stage):
    global quickprobs
    if killed_stage != 2:
        if killed_stage != 3:
            print("[MAIN STEP] Seperating Credible Regions...")
            reliable_regions = getReliableRegions(col_score, threshold, 0, 0, seq_file, dir_output)
            seperateReliableRegions(reliable_regions, real_output, dir_output)
            print("[MAIN STEP] Seperated Credible Regions.")
        else:
            killed_stage = 4
    else:
        killed_stage = 4
        os.system(quickprobs + "  " + seq_file + " > " + output_file)
    return killed_stage
