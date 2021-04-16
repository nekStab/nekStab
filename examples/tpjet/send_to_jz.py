from os import getcwd, system
from os.path import normpath, basename
import time

workdir = basename(normpath(getcwd()))
path = "rvpo014@jean-zay.idris.fr:/gpfswork/rech/vpo/rvpo014/2cyl/" + workdir + "/"
print(path)

file = []
file.append(['', '*'])

for i in file:
    folder = i.pop(0)
    print(folder)
    for j in i:
        print(j)
        system('rsync -auv '+path+folder+j+' .'+folder)
