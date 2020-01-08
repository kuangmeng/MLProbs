#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 13:19:32 2019

@author: mmkuang
"""
import os

def Refresh():
    os.system("rm -rf ./tmp/")
    os.system("mkdir ./tmp/")
    os.system("mkdir ./tmp/seperate_regions")

def Move(file_in, file_out):
    os.system("cp "+ file_in + " " + file_out)


