#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 00:17:17 2019

@author: mmkuang
"""
import os

def getTC_SP(qscore, output_file, reference_file):
    os.system(qscore + " -test " + output_file + " -ref " + reference_file + " > qscore.tc_sp")