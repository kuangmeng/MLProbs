# -*- coding: utf-8 -*-
"""
Spyder Editor

"""
import os
import time

def Preprocessing(bench, file, path):
    fileout = path +"/" + bench +"/"+ file
    filein = path +"/"+ bench +"/"+ file
    ret_list = []
    flag = 0
    filein_list = open(filein, 'r').read().splitlines()

    for line in filein_list:
        if len(line.strip()) > 0:
            if line.strip()[0] == '>':
                ret_list.append(line.strip())
                flag = 0
            else:
                if flag == 0:
                    ret_list.append(line.strip())
                    flag = 1
                else:
                    ret_list[len(ret_list) - 1] += line.strip()
    with open(fileout, 'w') as fileout_:
        for line in ret_list:
            fileout_.write(line + "\n")

def Run(prom, bench, file, path):
    Preprocessing(bench, file, path)
    f_t = os.path.exists("./output/" + bench)
    if not f_t:
        os.makedirs("./output/" + bench)
    os.system("%s %s/%s/%s ./output/%s/%s" % (prom, path, bench, file, bench, file))

def Compute(prom):
    path = "./BENCH"
    bench_list = os.listdir(path)
    tmp_bench = []
    total_time = []
    for bench in bench_list:
        if os.path.isdir(path + "/" + bench):
            tmp_bench.append(bench)
            tmp_time = 0.0
            file_num = 0
            file_list = os.listdir(path + "/" + bench)
            for fileidx in range(len(file_list)):
                tmp_str = ""
                if file_list[fileidx][0:1] == ".":
                    continue
                file_num += 1
                print("开始处理%d号文件: %s ..." % (fileidx, bench + "/" + file_list[fileidx]))
                start = time.time()
                Run(prom, bench, file_list[fileidx], path)
                file_write = open("train_write.txt", 'a')
                with open("./tmp/train_write.txt", 'r') as filein:
                    tmp_str = filein.read().splitlines()[0]

                end = time.time()
                tmp_time += (end - start)
                print("完成处理%d号文件！" % (fileidx))
                file_write.write(tmp_str + "\n")
                file_write.close()
            total_time.append(tmp_time / float(file_num))

    return tmp_bench, total_time

if __name__ == "__main__":
    tmp_bench, total_time = Compute("python MLProbs.py ")
    for i in range(len(tmp_bench)):
        print("%s: %f" % (tmp_bench[i], total_time[i]))
