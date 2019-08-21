#!/usr/bin/env python
import  xdrlib ,sys
import xlrd
from xlwt import *
import math

txt_file = "train_3.txt"
excel_file = "output.xlsx"

def open_excel(file):
    data = xlrd.open_workbook(file)
    return data

def ReadTXT(txt_file):
    ret_matrix = []
    with open(txt_file, 'r') as file_in:
        file_context = file_in.read().splitlines()
        for line in file_context:
            tmp_list = line.split("\t")
            ret_matrix.append(tmp_list)
    return ret_matrix

def ReadEXCEL(excel_file):
    ret_matrix = []
    data = open_excel(excel_file)
    table = data.sheets()[0]
    nrows = table.nrows
    best_ = 0.0
    best = [0 for i in range(2)]
    for j in range(nrows):
        row = table.row_values(j)

        if float(row[2]) >= float(row[1]):
            ret_matrix.append([row[0], 1])
            best[1] += 1
            best_ += float(row[1])
        else:
            ret_matrix.append([row[0], 0])
            best[0] += 1
            best_ += float(row[1])
    print(best)
    print(best_ / nrows)
    return ret_matrix

def getLabel(ret_matrix1, ret_matrix2):
    train = []
    label = []
    for i_line in ret_matrix1:
        for j_line in ret_matrix2:
            if i_line[0] == j_line[0]:
                train.append(i_line[1:])
                label.append(j_line[1])
                break
    return train, label

def PrepareData():

    ret_matrix1 = ReadTXT(txt_file)
    ret_matrix2 = ReadEXCEL(excel_file)
    train, label = getLabel(ret_matrix1, ret_matrix2)

    para = [0.0 for i in range(len(train[0]) * 2)]
    train_final = []
    for i in range(len(train[0])):
        for j in range(len(train)):
            if float(train[j][i]) > para[i * 2]:
                para[i * 2] = float(train[j][i])
            if float(train[j][i]) < para[i * 2 + 1]:
                para[i * 2 + 1] = float(train[j][i])

    for i in range(len(train)):
        tmp_arr = []
        for j in range(len(train[0])):
            tmp_arr.append((float(train[i][j]) - para[j * 2 + 1]) / (para[j * 2] - para[j * 2 + 1]))
        train_final.append(tmp_arr)

    with open("para.txt", 'w') as fileout:
        for i in para:
            fileout.write(str(i) + "\n")

    return train_final, label
