from os import getcwd, system
from os.path import normpath, basename
import time

workdir = basename(normpath(getcwd()))
path = "rfrantz@minerve:/home/rfrantz/tpjet_report/" + workdir + "/"
print(path)

file = []
file.append(['', '*'])
file.append(['/figures', '*'])
file.append(['', 'Spectre_*','*.his','*.py', 'Re*', 'KRY*', 'PER*'])

for i in file:
    folder = i.pop(0)
    print(folder)
    for j in i:
        print(j)
        system('rsync -auv '+path+folder+j+' .'+folder)
