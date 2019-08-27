#!/usr/bin/env python

import os
from xlwt import *
import sys

def ReadList(dir):
    tmp_bench_list = os.listdir("./result")
    bench_list = []
    for item in tmp_bench_list:
        if not item[0] == '.':
            bench_list.append(item)
    getTCmatrix(dir, bench_list)

def getTCmatrix(dir, bench_list):
    TCmatrix = []
    for bench in bench_list:
        file_list = os.listdir("./result/%s"%(bench))
        for file in file_list:
            if not file[0] == '.':
                row = []
                row.append(bench + "." + file)
                files = "./result/" + bench + "/" + file
                if not os.path.exists(files):
                    continue
                elif not os.path.getsize(files):
                    continue
                filein = open(files)
                file_context = filein.read().splitlines()
                filein.close()
                if len(file_context) == 0:
                    tmp_q = tmp_tc = 0
                else:
                    tmp_list = file_context[0].split(";")
                    tmp_tc = float(tmp_list[3].split("=")[1])
                    tmp_q = float(tmp_list[2].split("=")[1])
                row.append(tmp_tc)
                row.append(tmp_q)
            TCmatrix.append(row)
        toExcel(TCmatrix, dir, bench)
        TCmatrix = []

def toExcel(matrix, dir, name):
    file = Workbook(encoding = 'utf-8')
    #指定file以utf-8的格式打开
    table = file.add_sheet('data')
    style = XFStyle()
    style.num_format_str = '0.000000'
    table.write(0, 1, "TC-score")
    table.write(0, 2, "SP-score")
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            table.write(i + 1, j, matrix[i][j], style)
    file.save("./" + dir + "/" + name + '.xlsx')

def Run(bench, file, dir):
    files = "./result/%s"%(bench)
    if not os.path.exists(files):
        os.mkdir(files)
    Process("./ref/" + bench + "/" + file)
    Process("./" + dir + "/" + bench + "/" + file)
    os.system("./qscore/qscore -test ./%s/%s/%s -ref ./ref/%s/%s > ./result/%s/%s"%(dir, bench, file, bench, file, bench, file))

def Process(file):
    filein = open(file, 'r')
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
    lens = len(value)
    dickeys = sorted(dic.keys())
    with open(file, 'w') as fileout:
        for key_ in dickeys:
            fileout.write(key_ + "\n")
            fileout.write(dic[key_] + "\n")

def Compute(dir):
    bench_list = os.listdir("./" + dir)
    for item in bench_list:
        if not item[0] == '.':
            file_list = os.listdir("./" + dir + "/" + item)
            for file in file_list:
                if not file[0] == '.':
                    Run(item, file, dir)

if __name__ == "__main__":
    dir = sys.argv[1]
    Compute(dir)
    ReadList(dir)
