import os

bench = "../BENCH/"

quick = "./quickprobs "

output = "./output/"

def getOutput_Dir(bench):
    bench_list = os.listdir(bench)
    for bench_ in bench_list:
        if os.path.isdir(bench + bench_):
            getOutput_File(bench, bench_)

def getOutput_File(bench, bench_):
    if not os.path.exists(output + bench_):
        os.mkdir(output + bench_)
    file_list = os.listdir(bench + bench_)
    for file in file_list:
        if file[0] != '.':
            os.system(quick + bench + bench_ + "/" + file + " > " + output + bench_ + "/" + file)

if __name__ == '__main__':
    getOutput_Dir(bench)
