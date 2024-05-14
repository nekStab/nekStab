#bin/python3

from statistics import geometric_mean
import numpy as np
import subprocess, sys, os
sys.path.append("..")
from autotools import *

curv_radius = 10.0 / 3.0
print('curv_radius = ', curv_radius)
slices_per_length = np.ceil(100 / (2.0 * np.pi * curv_radius))
print('slices_per_length = ', slices_per_length)
length = 1.0
print('length = ', length)
nslices = int (length * slices_per_length)
print('nslices = ', nslices)
nelf = 108 
print('nelf = ', nelf)
dp0ds = 0.0794600
print('dp0ds = ', dp0ds)
dp1dsr = 0.0
print('dp1dsr = ', dp1dsr)
re = 1700.0
print('re = ', re)
lx1 = 6  # 6, 8, 10
print('lx1 = ', lx1)
lelg = nelf * nslices
print('lelg = ', lelg)
case = 'torus'
print('case = ', case)

# read the file into a list of lines
with open('geom/templates/n2to3_cmd_template.txt', 'r') as file:
    lines = file.readlines()
lines[3] = str(nslices) + '\n'  # index is 3 because we start counting from 0
with open('geom/templates/n2to3_cmd_template.txt', 'w') as file:
    file.writelines(lines)

if nelf == 108:
    # os.chdir('geom')
    # subprocess.call(["./mkmsh.sh", "torus_coarse2D"])
    # os.chdir('..')
    subprocess.call(["./set_mesh.sh", "torus_coarse3D"])
elif nelf == 225:
    subprocess.call(["./set_mesh.sh", "torus_medium3D"])
elif nelf == 305:
    subprocess.call(["./set_mesh.sh", "torus_fine3D"])

# read the file into a list of lines
with open('./HELIXD', 'r') as file:
    lines = file.readlines()

# replace the lines with the new values
lines[0] = f"      real, parameter :: curv_radius = {curv_radius}d0\n"
lines[1] = f"      real, parameter :: length = {length}d0\n"
lines[2] = f"      integer, parameter :: nslices = {nslices}\n"
lines[3] = f"      integer, parameter :: nelf = {nelf}\n"
lines[4] = f"      real, parameter :: dp0ds = {dp0ds}d0\n"
lines[5] = f"      real, parameter :: dp1dsr = {dp1dsr}d0\n"

# write the lines back to the file
with open('./HELIXD', 'w') as file:
    file.writelines(lines)

lxd = int(lx1 * 1.5)
    #   parameter (lx1=06)!06,08,10,12   ! GLL points per element along each direction
    #   parameter (lxd=09)!09,12,15,18   ! GL  points for 

# read the file into a list of lines
with open('./SIZE', 'r') as file:
    lines = file.readlines()

# replace lines 14 and 15 with the new lines
lines[12] = f"      parameter (lx1={lx1:02}) ! GLL points per element along each direction\n"
lines[13] = f"      parameter (lxd={lxd:02}) ! GL  points\n"
lines[14] = f"      parameter (lelg={lelg}) ! max number of global elements \n"





# write the lines back to the file
with open('./SIZE', 'w') as file:
    file.writelines(lines)

params = [
    # ('GENERAL', 'startFrom', bf),
    # ('GENERAL', 'endTime', '1.25'),
    # ('GENERAL', 'userParam01', '2'),
    # ('GENERAL', 'userParam07', '200'),
    # ('GENERAL', 'userParam10', '1.7'),
    # ('PRESSURE', 'residualtol', str(tol)),
    # ('VELOCITY', 'residualtol', str(tol)),
    ('VELOCITY', 'viscosity', f"-{float(re)}")
]
for section, key, value in params:
    c_pf(case+'.par', case+'.par', {section: {key: value}})
