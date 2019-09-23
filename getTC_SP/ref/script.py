import os 

filelist = []

with open('testlist.txt', 'r') as filein:
    filelist = filein.read().splitlines()

for file in filelist:
    os.system("cp %s %s"%(file, "new_test/" + file.split('/')[1]))