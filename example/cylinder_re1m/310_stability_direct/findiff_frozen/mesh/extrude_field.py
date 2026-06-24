#!/usr/bin/env python
#
#
# jcanton@gmail.com

import numpy as np
# install with `pip install (--upgrade) pymech`
from pymech import neksuite as ns

q2dr = ns.readnek('190_2d/BF_1cyl0.f00001') # saturated 2d solution
field = q2dr
print(field.endian, field.istep, field.lr1, field.nbc, field.ncurv, field.ndim, field.nel, field.time, field.var, field.wdsz,field)

# initialize extruded field
q3dr = ns.readnek('190_sym/1cyl0.f00001') # 3d mesh
field = q3dr
#print(field.endian, field.istep, field.lr1, field.nbc, field.ncurv, field.ndim, field.nel, field.time, field.var, field.wdsz,field)

q3dr.time   = q2dr.time
q3dr.istep  = 0

#------------------------------------------------------------------------------
# extrude
#
for iel in range(q3dr.nel):

    if (iel%q2dr.nel == 0):
        print('\telem {0:d}'.format(iel))

    for iz in range(q3dr.lr1[2]):
        for iy in range(q3dr.lr1[1]):
            for ix in range(q3dr.lr1[0]):

                # # coordinates from 2D field (if needed, probably not)
                # x = q2dr.elem[iel % q2dr.nel].pos[0, 0, iy, ix]
                # y = q2dr.elem[iel % q2dr.nel].pos[1, 0, iy, ix]
                # z = q2dr.elem[iel % q2dr.nel].pos[2, 0, iy, ix]
                # #
                # q3dr.elem[iel].pos[0, iz, iy, ix] = x
                # q3dr.elem[iel].pos[1, iz, iy, ix] = y
                # q3dr.elem[iel].pos[2, iz, iy, ix] = z

                # velocity from 2D field
                vx = q2dr.elem[iel % q2dr.nel].vel[0, 0, iy, ix]
                vy = q2dr.elem[iel % q2dr.nel].vel[1, 0, iy, ix]
                vz = 0.000 #q2dr.elem[iel % q2dr.nel].vel[2, 0, iy, ix]
                p  = q2dr.elem[iel % q2dr.nel].pres[0,0, iy, ix]
                t  = q2dr.elem[iel % q2dr.nel].temp[0,0, iy, ix]
                #
                q3dr.elem[iel].vel[0, iz, iy, ix] = vx
                q3dr.elem[iel].vel[1, iz, iy, ix] = vy
                q3dr.elem[iel].vel[2, iz, iy, ix] = vz
                q3dr.elem[iel].pres[0, iz, iy, ix] = p
                q3dr.elem[iel].temp[0, iz, iy, ix] = t

#little 115592 [6, 6, 1] 0 [] 2 1480 1000.0 [2, 2, 1, 1, 0] 4
#little 0 [6, 6, 6] 0 [] 3 14800 0.0 [3, 3, 1, 1, 0] 4


#------------------------------------------------------------------------------
# save
#
filename = '190_sym/BFEXT_1cyl0.f00001'
print('Saving extuded solution to',filename)
ns.writenek(filename, q3dr)
