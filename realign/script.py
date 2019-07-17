import os
import time

bench = "../bench_all/"

quick = "./quickprobs "

output = "./output/"

def getOutput_Dir(bench):
    bench_list = os.listdir(bench)
    fileout = open("time.txt", 'w')
    for bench_ in bench_list:
        if os.path.isdir(bench + bench_):
            start = time.time()
            lens = getOutput_File(bench, bench_)
            end = time.time()
            times = (end - start) / lens
            fileout.write(bench_ + str(times) + "\n")


def getOutput_File(bench, bench_):
    if not os.path.exists(output + bench_):
        os.mkdir(output + bench_)
    file_list = os.listdir(bench + bench_ + "/in/")
    lens = 0
    for file in file_list:
        if file[0] != '.':
            lens += 1
            os.system(quick + bench + bench_ + "/in/" + file + " > " + output + bench_ + "/" + file)
            print("Processing %d"%(lens))
    return lens

if __name__ == '__main__':
    getOutput_Dir(bench)
