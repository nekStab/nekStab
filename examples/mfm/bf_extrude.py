# --> Import standard python libraries.
import numpy as np

import struct
import matplotlib.pyplot as plt

params = {'text.usetex' : False,
          'font.size' : 8,
		  }
plt.rcParams['figure.facecolor'] = '1'
plt.rcParams.update(params)

plt.style.use('seaborn-white')

fig_width = 5.33

# --> Define the plotting function.
def discrete_cmap(N=9, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    if N == 9:
        color_list[4] = (1., 1., 1., 1.)
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

cmap = discrete_cmap(9, base_cmap=plt.cm.hot)

def readnek(fname):

	"""
	Utility function to read nek5000 binary files.

	Parameters
	----------
	fname : string
    	Name of the file to be read.

	Returns
	-------
	type
    	Description of returned object.

	"""

	####################################
	#####     Parse the header     #####
	####################################

	# --> Open the file.
	try:
	    infile = open(fname, "rb")
	except IOError as e:
	    print("I/O erro({0}): {1}".format(e.errno, e.strerror))
	    return -1

	# --> Read header.
	header = infile.read(132).split()

	# --> Get word size.
	wdsz = int(header[1])
	if wdsz == 4:
	    realtype = "f"
	elif wdsz == 8:
	    readltype = "d"
	else:
	    print("ERROR: Could not interpret real type (wdsz = %i)" %wdsz)
	    return -2

	# --> Get polynomial order.
	lr1 = [int(i) for i in header[2:5]]

	# --> Compute total number of points per element.
	npel = np.prod(lr1)

	# --> Number of physical dimensions.
	ndim = 2 + (lr1[2]>1)

	# --> Get the number of elements.
	nel = int(header[5])

	# --> Get the number of elements in that file.
	nelf = int(header[6])

	# --> Get current time.
	time = float(header[7])

	# --> Get current time-step.
	istep = int(header[8])

	# --> Get file ID.
	fid = int(header[9])

	# --> Get total number of files.
	nf = int(header[10])

	# --> Get variables [XUPT]
	vars = header[11].decode("utf-8")
	var = np.zeros(5, dtype=np.int)
	for v in vars:
	    if v == "X":
	        var[0] = ndim
	    elif v == "U":
	        var[1] = ndim
	    elif v == "P":
	        var[2] = 1
	    elif v == "T":
	        var[3] = 1
	    elif v == "S":
	    	v[4] = 0

	# --> Total number of scalar fields to be read.
	nfields = var.sum()

	# --> Identify endian encoding.
	etagb = infile.read(4)
	etagL = struct.unpack('<f', etagb)[0]
	etagL = int(etagL*1e5)/1e5
	etagB = struct.unpack('>f', etagb)[0]
	etagB = int(etagB*1e5)/1e5

	if etagL == 6.54321:
		emode = '<'
	elif etagB == 6.54321:
		emode = '>'
	else:
		print('ERROR: could not interpret endianness')
		return -3

	# --> Read the element map.
	elmap = infile.read(4*nelf)
	elmap = list(struct.unpack(emode+nelf*"i", elmap))

	########################################
	#####     Read the actual data     #####
	########################################

	# --> Numpy array container for the data.
	data = np.zeros((nelf, npel, nfields))

	for ivar in range(len(var)):
		idim0 = sum(var[:ivar])
		for iel in elmap:
			for idim in range(var[ivar]):
				x = infile.read(npel*wdsz)
				data[iel-1, :, idim+idim0] = struct.unpack(emode+npel*realtype, x)

	# --> Close the file.
	infile.close()

	# --> Output dictionnary.
	output = {
		"data": data,
		"lr1": lr1,
		"elmap": elmap,
		"time": time,
		"istep": istep,
		"emode": emode,
		"wdsz": wdsz,
		"header": header,
		"fields": vars
	}

	return output

def writenek(fname, dico):

	##################################
	#####     Initialization     #####
	##################################

	# --> Open the file.
	try:
		outfile = open(fname, "wb")
	except IOError as e:
		print("I/O erro ({0}) : {1}".format(e.errno, e.strerror))
		return -1

	# --> Misc.
	fid, nf = 0, 1

	# --> Polynomial order.
	lr1 = dico["lr1"]

	# --> Get data.
	data = dico["data"]

	# --> Get element map.
	elmap = dico["elmap"]

	# --> Number of elements, points per element.
	nel, npel, _ = data.shape
	assert npel == np.prod(lr1)

	# --> Number of active dimensions.
	ndim = 2 + (lr1[2]>1)

	# --> Get fields to be written.
	var = np.zeros(5, dtype=np.int)
	for v in dico["fields"]:
	    if v == "X":
	        var[0] = ndim
	    elif v == "U":
	        var[1] = ndim
	    elif v == "P":
	        var[2] = 1
	    elif v == "T":
	        var[3] = 1
	    elif v == "S":
	    	v[4] = 0

	# --> Get word size.
	if dico["wdsz"] == 4:
		realtype = "f"
	elif dico["wdsz"] == 8:
		realtype = "d"
	else:
	    print("ERROR: Could not interpret real type (wdsz = %i)" %wdsz)
	    return -2

	# --> Get endianness.
	emode = dico["emode"]

	# --> Generate header.
	header = "#std %1i %2i %2i %2i %10i %10i %20.13E %9i %6i %6i %s\n" %(
		dico["wdsz"], lr1[0], lr1[1], lr1[2], nel, nel,
		dico["time"], dico["istep"], fid, nf, dico["fields"])

	header = header.ljust(132).encode("utf-8")

	# --> Write header.
	outfile.write(header)

	etagb = struct.pack(emode+"f", 6.54321)
	outfile.write(etagb)

	outfile.write(struct.pack(emode+nel*"i", *elmap))

	##############################
	#####     Write data     #####
	##############################

	for ivar in range(len(var)):
		idim0 = np.sum(var[:ivar])
		for iel in elmap:
			for idim in range(var[ivar]):
				x = struct.pack(emode+npel*realtype, *data[iel-1, :, idim+idim0])
				outfile.write(x)

	# --> Close the file.
	outfile.close()
	return

def pltmd(x, y, q, savename=None, qual=300):

    # --> Create the triangulation.
    from matplotlib.tri import Triangulation
    triang = Triangulation(x, y)

    # --> Mask the cylinder.
    xmid = x[triang.triangles].mean(axis=1)
    ymid = y[triang.triangles].mean(axis=1)

    mask = np.where(xmid**2 + ymid**2 < 0.5**2, 1, 0)
    triang.set_mask(mask)

    # --> Plot the figure.
    fig = plt.figure(figsize=(fig_width/2, fig_width/4))
    ax = fig.gca()

    # --> Shaded contour plot of the POD mode.
    max_val = abs(q).max()
    ax.tricontourf(triang, q,cmap=plt.cm.hot)

    # --> Line contour.
    levels = np.array([0.2,0.4, 0.5, 0.75])*max_val
    ax.tricontour(triang, q, levels=levels, colors='k', linewidths=0.5)
    ax.tricontour(triang, q, levels=-np.flipud(levels), colors='k', linewidths=0.5)

    # --> Axes related.
    ax.set_yticks([])
    ax.set_ylim(-5., 5.)
    ax.set_xlim(-2.5, 47.5)
    ax.set_xlabel(r'$x$')
    ax.set_aspect('equal')
    ax.locator_params(axis='x', nbins=4.)

    # --> Hide the cylinder.
    from matplotlib.patches import Circle
    patch = Circle(xy=(0, 0), radius=0.5, facecolor='lightgray', edgecolor='black', lw=0.5)
    ax.add_patch(patch)

    if savename is not None:
        plt.savefig(savename, bbox_inches='tight', dpi=qual)

    plt.close()
    plt.clf()
    return


if __name__ == "__main__":

    q2dr = readnek('BF_mfm0.f00001') # 2d solution
    q2d = q2dr["data"]
    print(q2d.shape)

    q3dr = readnek('mfm3d0.f00001') # 3d solution at 1 iteration
    q3d = q3dr["data"]
    print(q3d.shape)

    lr1 = q2dr["lr1"][0]

    for z in range(nz):
    	q3d[:,:,z,2] = q2d[:,:,2] #u from 2D
    	q3d[:,:,z,3] = q2d[:,:,3] #v from 2D
    	q3d[:,:,z,4] = 0.         #force zero
    	q3d[:,:,z,5] = q2d[:,:,4] #p from 2D
    	q3d[:,:,z,6] = q2d[:,:,5] #p from 2D

    writenek('BF_mfm3d0.f00001', q3d) # save 3d extruded solution
