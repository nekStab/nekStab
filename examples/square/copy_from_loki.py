from os import getcwd, system
from os.path import normpath, basename
import time

workdir = basename(normpath(getcwd()))
path = "acastriotta@loki:/home/acastriotta/Re50/"
print(path)

file = []

file.append(['', 'Square0.f02562'])
file.append(['', 'Square.his'])

for i in file:
    folder = i.pop(0)
    print(folder)
    for j in i:
        print(j)
        system('rsync -auv '+path+folder+j+' .'+folder)
