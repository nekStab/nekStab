
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
import numpy as np
from scipy import interpolate
import struct

params = {'text.usetex' : False,
          'font.size' : 8,
          #'font.family' : 'lmodern'
          }
#plt.rcParams['text.latex.preamble']=[r'\usepackage{lmodern}']
plt.rcParams['figure.facecolor'] = '1'
plt.rcParams.update(params)
plt.style.use('seaborn-white')

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""
    from matplotlib.colors import ListedColormap
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return ListedColormap(color_list, name=cmap_name)

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

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

def plt2d(x,y,u,v,   q, savename=None):

    # --> Create the triangulation.
    from matplotlib.tri import Triangulation
    triang = Triangulation(x, y)

    # --> Mask the cylinder.
    if ifobjmask:
        xmid = x[triang.triangles].mean(axis=1)
        ymid = y[triang.triangles].mean(axis=1)
        mask = np.where(xmid**2 + ymid**2 < 0.5**2, 1, 0)
        triang.set_mask(mask)

    # --> Plot the figure.
    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.gca()

    # --> Shaded contour plot of the mode.
    mxv = abs(q).max()
    mnv = abs(q).min()

    if normalize:
        #q /= mxv
        mxv = abs(q).max()
        mnv = abs(q).min()
        tcf = ax.tripcolor(triang, q, shading='gouraud', cmap=cmap, norm=MidpointNormalize(midpoint=0.000))
        if colorbar:
            fig.colorbar(tcf,extend='neither', orientation='vertical',ticks=[-1, 0, 1])
            cbar.ax.set_yticklabels(['-1', '0', '1'])  # vertically oriented colorbar

        #neither' | 'both' | 'min' | 'max'
    else:
        tcf = ax.tricontourf(triang, q, cmap=cmap)
        if colorbar:
            fig.colorbar(tcf,extend='neither', orientation='vertical')
            #neither' | 'both' | 'min' | 'max'


    # --> Line contour.
    if ifiso:
        levels = np.linspace(iso_in*mxv,iso_fi*mxv,iso_st)
        #levels = np.array([0.1,0.2,0.5,0.8,0.9])*max_val
        ax.tricontour(triang, q, levels=levels, colors=isocolor, linewidths=isoline)
        ax.tricontour(triang, q, levels=-np.flipud(levels), colors=isocolor, linewidths=isoline)

    if ifstream:
        nx = int((xmx-xmn)*mult)
        ny = int((ymx-ymn)*mult)
        x2 = np.linspace(xmn,xmx, nx, endpoint=True)
        y2 = np.linspace(ymn,ymx, ny, endpoint=True)
        x2, y2 = np.meshgrid(x2, y2)
        trg = np.c_[x2.ravel(),y2.ravel()]
        src = np.c_[x,y]
        u2 = interpolate.NearestNDInterpolator(src, u)(trg)
        v2 = interpolate.NearestNDInterpolator(src, v)(trg)
        u2 = np.reshape(u2, (ny, nx))
        v2 = np.reshape(v2, (ny, nx))
        magnitude = np.sqrt(u2**2 + v2**2)
        ax.streamplot(x2, y2, u2, v2, color=magnitude, cmap=plt.cm.hot, density=density,linewidth=linewidth, \
        arrowsize=arrowsize, minlength=minlength, maxlength=maxlength)

    # --> Axes related.
    ax.set_yticks([])
    if 'xmx' in globals():
        ax.set_ylim(ymn,ymx)
        ax.set_xlim(xmn,xmx)
    else:
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$r$')
    #ax.set_aspect('equal')
    ax.locator_params(axis='x', nbins=4.)

    # --> Hide the cylinder.
    if ifobjmask:
        from matplotlib.patches import Circle
        patch = Circle(xy=(0, 0), radius=0.5, facecolor='lightgray', edgecolor='black', lw=0.5)
        ax.add_patch(patch)

    if savename is not None:
        plt.savefig(savename, bbox_inches='tight', dpi=qual)
    plt.close();plt.clf()
    return

def plotm(ft,start,end,incr=1,*args):
    for i in range(start, end, incr):
        fname = ft + case + '0.f' + str(i).zfill(5)
        #try:
        print('file=',fname)
        q = readnek(fname)["data"]
        x = q[:,:,0].ravel()
        y = q[:,:,1].ravel()

        if if3d:
            u = q[:,:,4].ravel(); v = q[:,:,5].ravel()
            z = q[:,:,3].ravel(); w = q[:,:,6].ravel()
            m = np.sqrt(u**2 + v**2 + w**2)
        else:
            u = q[:,:,2].ravel()
            v = q[:,:,3].ravel()
            m = np.sqrt(u**2 + v**2)

        plt2d(x,y,u,v, u, savename=outfldr+'/'+ft+'/'+'u_'+ft+'_'+case+'_'+str(i).zfill(5))
        plt2d(x,y,u,v, v, savename=outfldr+'/'+ft+'/'+'v_'+ft+'_'+case+'_'+str(i).zfill(5))
        plt2d(x,y,u,v, m, savename=outfldr+'/'+ft+'/'+'M_'+ft+'_'+case+'_'+str(i).zfill(5))

        #except:
        #    print('file=',fname,'not found!')
        #    pass

if __name__ == "__main__":

    case = 'tpjet'
    fldr = './'
    outfldr = fldr + '/figures/'
    #cmap = discrete_cmap(9, 'RdBu_r')
    #'Perceptually Uniform Sequential'] = ['viridis', 'plasma', 'inferno', 'magma', 'cividis']
    #'Sequential'] = ['Greys', ']
    #'Sequential (2)'] = ['binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink', 'summer', 'autumn', 'winter', 'cool', 'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper']
    #'Diverging'] = ['Spectral', 'coolwarm', 'bwr', 'seismic','PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu','RdYlBu', 'RdYlGn',]
    #'Cyclic'] = ['twilight', 'twilight_shifted', 'hsv']

    if3d = False
    ifobjmask = False
    fig_width = 5.33
    fig_height = fig_width/7.5
    qual = 400

    xmn = 0.0
    xmx = 40.0
    ymn = 0.0
    ymx = 2.0
    isoline = 0.5
    isocolor = 'k'
    iso_in = 0.10
    iso_fi = 0.90
    iso_st = 4

    mult = 20
    density = 0.9
    linewidth = 0.7
    arrowsize = 0.8
    minlength = 1
    maxlength = 10

    cmap = discrete_cmap(23,'RdBu') #peturbations with diverging colormaps
    normalize = True
    colorbar = False
    ifiso = True
    ifstream = False
    #plotm('Re_',1,3)
    #plotm('Rev',1,3)
    plotm('PER',1,10)
    plotm('',,10)
