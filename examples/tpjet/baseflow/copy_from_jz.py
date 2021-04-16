from os import getcwd, system
from os.path import normpath, basename
import time

workdir = basename(normpath(getcwd()))
path = "rvpo014@jean-zay.idris.fr:/gpfswork/rech/vpo/rvpo014/tpjet/tpjet_2000_tdf/"# + workdir + "/"
print(path)

file = []

#file.append(['', 'SIZE', 'x_*','*.py','*.sh','makefile_usr.inc'])
#file.append(['', '*.usr','*.par'])
#file.append(['', 'tpjet0.f02401','Revtpjet0.f0000*'])
file.append(['', 'tpjet.box'])
#file.append(['', 'Re_tpjet0.f0000*','Im_tpjet0.f0000*','PERtpjet0.f0000*'])
#file.append(['', 'Spectre*'])
#file.append(['', 'HEStpjet*'])


for i in file:
    folder = i.pop(0)
    print(folder)
    for j in i:
        print(j)
        system('rsync -auv '+path+folder+j+' .'+folder)
