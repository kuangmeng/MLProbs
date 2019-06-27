# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import sys

def getLens(filename):
    dic = {}
    has_key = False
    value = ""
    key = ""
    filein = open(filename, 'r')
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
    
    for item in dic.keys():
        print("%s: %d" % (item, len(dic[item])))

if __name__ == "__main__":
    getLens(sys.argv[1])