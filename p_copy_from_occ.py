from os import getcwd, system
from os.path import normpath, basename
import time

workdir = basename(normpath(getcwd()))
path = "rfrantz@occigen.cines.fr:/scratch/cnt0027/ldf6362/rfrantz/" + workdir + "/"
print(path)

print("ck4k3yzkvhr2")

file = []
file.append(['','*'])

for i in file:
    folder = i.pop(0)
    print(folder)
    for j in i:
        print(j)
        system('rsync -auv '+path+folder+j+' .'+folder)
