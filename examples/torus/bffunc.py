#!/usr/bin/env python
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.tri import Triangulation
from scipy.spatial import Delaunay
from scipy.ndimage import gaussian_filter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from django.urls import path, include
import struct, os, timeit

# bhm
# classic
# grayscale
# seaborn-paper
# seaborn-colorblind
# plt.style.use("seaborn-v0_8-colorblind")
# plt.style.use("seaborn-v0_8-ticks")

#conda install scienceplots
#pip install SciencePlots
import scienceplots 
plt.style.use("science")

# http://tonysyu.github.io/raw_content/matplotlib-style-gallery/gallery.html
'''
RcParams({'axes.axisbelow': True,
          'axes.edgecolor': 'white',
          'axes.facecolor': '#EAEAF2',
          'axes.grid': True,
          'axes.labelcolor': '.15',
          'axes.linewidth': 0.0,
          'figure.facecolor': 'white',
          'font.family': ['sans-serif'],
          'font.sans-serif': ['Arial',
                              'Liberation Sans',
                              'DejaVu Sans',
                              'Bitstream Vera Sans',
                              'sans-serif'],
          'grid.color': 'white',
          'grid.linestyle': '-',
          'image.cmap': 'Greys',
          'legend.frameon': False,
          'legend.numpoints': 1,
          'legend.scatterpoints': 1,
          'lines.solid_capstyle': 'round',
          'text.color': '.15',
          'xtick.color': '.15',
          'xtick.direction': 'out',
          'xtick.major.size': 0.0,
          'xtick.minor.size': 0.0,
          'ytick.color': '.15',
          'ytick.direction': 'out',
          'ytick.major.size': 0.0,
          'ytick.minor.size': 0.0})
'''

#https://www.cambridge.org/core/services/authors/journals/journals-artwork-guide
plt.rcParams.update({'text.usetex': False,
                     'font.size': 8, # for JFM 
                     'legend.fontsize': 8,
                     'legend.handlelength': 1.})

formt = 'png' # easy to change in final version
ajust = 'tight'
qual = 900 # DPI
symbols = ['o', 'v', '^', 's', 'd', 'x', '*', '+', '.']
colors = ['royalblue', 'seagreen', 'orange', 'yellow', 'm', 'brown', 'gray' ] # old colors 
colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00'] # colorblind friendly
# colors = {
#     'blue':   [55,  126, 184],  #377eb8 
#     'orange': [255, 127, 0],    #ff7f00
#     'green':  [77,  175, 74],   #4daf4a
#     'pink':   [247, 129, 191],  #f781bf
#     'brown':  [166, 86,  40],   #a65628
#     'purple': [152, 78,  163],  #984ea3
#     'gray':   [153, 153, 153],  #999999
#     'red':    [228, 26,  28],   #e41a1c
#     'yellow': [222, 222, 0]     #dede00
# }; opacity = 0.8
# c_str = {k:f'rgba({v[0]},{v[1]},{v[2]},{opacity})'
#          for (k, v) in colors.items()}
# c_str['yellow'] # Gives the rgba string for 'yellow'

## \textheight = 674pt
## \textwidth = 426pt
## single column 240.7103 pt

# fig_width = 426/72.27
# fig_height = (16/9)*fig_width # nice aspect ratio

fig_width = 6.32
fig_height = fig_width

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""
    from matplotlib.colors import ListedColormap
    base = mpl.colormaps.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return ListedColormap(color_list, name=cmap_name)

# class MidpointNormalize(colors.Normalize):
#     def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
#         self.midpoint = midpoint
#         colors.Normalize.__init__(self, vmin, vmax, clip)

#     def __call__(self, value, clip=None):
#         x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
#         return np.ma.masked_array(np.interp(value, x, y))

class Spectre(object):
    def __init__(self, filename):
        #data = np.transpose(np.loadtxt(filename))
        data = np.transpose(np.genfromtxt(filename))
        self.vp_real  = data[0]
        self.vp_imag  = data[1]
        self.residu   = data[2]
        del data

class Spectre_c(object):
    def __init__(self, filename):
        #data = np.transpose(np.loadtxt(filename))
        data = np.transpose(np.genfromtxt(filename))
        self.vp_real  = data[0]
        self.vp_imag  = data[1]
        del data

def readnek(fname):
    
    try:
        infile = open(fname, "rb")
    except IOError as e:
        print("I/O erro({0}): {1}".format(e.errno, e.strerror))
        return -1

    header = infile.read(132).split()

    # --> Get word size.
    wdsz = int(header[1])
    if wdsz == 4:
        realtype = "f"
    elif wdsz == 8:
        realtype = "d"
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
    var = np.zeros(5, dtype=int)
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
            var[4] = 1

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

def plot_BF_cube_sphere(savename,qc,qs,lim):
    
    xmin = -2.8; xmax = 9.
    ymin = -2; ymax = 2
    zz=64
    cmap = discrete_cmap(zz,'gray')
    
    fig,ax = plt.subplots(ncols=2, nrows=2)
    fig.set_size_inches(fig_width,fig_height/2.222)
    
    ax[0,0].axis('equal')
    ax[0,1].axis('equal')
    ax[1,0].axis('equal')
    ax[1,1].axis('equal')
      
    ax[0,0].set_ylim(ymin,ymax)
    ax[0,1].set_ylim(ymin,ymax)
    ax[1,0].set_ylim(ymin,ymax)
    ax[1,1].set_ylim(ymin,ymax)
        
    ax[0,0].set_xlim(xmin,xmax) 
    ax[0,1].set_xlim(xmin,xmax)
    ax[1,0].set_xlim(xmin,xmax) 
    ax[1,1].set_xlim(xmin,xmax)
    
    ax[0,0].set_ylabel(r'$y$',labelpad=-0) 
    ax[1,0].set_ylabel(r'$z$',labelpad=-0)
    ax[1,0].set_xlabel(r'$x$',labelpad=-0)
    ax[1,1].set_xlabel(r'$x$',labelpad=-0)
    
    ax[0,0].set_xticklabels([])
    ax[0,1].set_yticklabels([])
    ax[0,1].set_xticklabels([])
    ax[1,1].set_yticklabels([])

    levels = np.linspace(lim[0], lim[1], zz)
    lvs = [0,1] # np.linspace(lim[0]*0.8, lim[1]*0.8,4)

    ##### PLOT A #####
    axp0 = ax[0,0].tricontourf(qc['y']['tr'], qc['y']['u'], cmap=cmap,levels=levels,extend='both')
    ax[0,0].tricontour(qc['y']['tr'], qc['y']['u'], levels=lvs,  colors='k', linewidths=0.2)
    ax[0,0].tricontour(qc['y']['tr'], qc['y']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes0 = inset_axes(ax[0,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$u$', fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)

    axp1 = ax[1,0].tricontourf(qc['z']['tr'], qc['z']['u'], cmap=cmap,levels=levels,extend='both')
    ax[1,0].tricontour(qc['z']['tr'], qc['z']['u'], levels=lvs,  colors='k', linewidths=0.2)
    ax[1,0].tricontour(qc['z']['tr'], qc['z']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes1 = inset_axes(ax[1,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$u$', fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
    
    ##### PLOT B #####
    axp2 = ax[0,1].tricontourf(qs['y']['tr'], qs['y']['u'], cmap=cmap,levels=levels,extend='both')
    ax[0,1].tricontour(qs['y']['tr'], qs['y']['u'], levels=lvs,  colors='k', linewidths=0.2)
    ax[0,1].tricontour(qs['y']['tr'], qs['y']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes2 = inset_axes(ax[0,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp2,orientation='horizontal',cax=cbaxes2)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$u$', fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
    
    axp3 = ax[1,1].tricontourf(qs['z']['tr'], qs['z']['u'], cmap=cmap,levels=levels,extend='both')
    ax[1,1].tricontour(qs['z']['tr'], qs['z']['u'], levels=lvs,  colors='k', linewidths=0.2)
    ax[1,1].tricontour(qs['z']['tr'], qs['z']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes3 = inset_axes(ax[1,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp3,orientation='horizontal',cax=cbaxes3)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$u$', fontsize=7,labelpad=-8)
    cbar.set_ticklabels(ticks)
    
    ax[0,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[0,1].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,1].grid(True, c='gray', ls=':', lw=0.6)
    
    ##### Patches #####
    pc1 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    pc2 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    ps1 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    ps2 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    
    ax[0,0].add_artist(pc1)
    ax[1,0].add_artist(pc2)
    ax[0,1].add_artist(ps1)
    ax[1,1].add_artist(ps2)
    
    # Move patches to the front
    pc1.set_zorder(2)
    ps1.set_zorder(2)
    
    ax[0,0].set_title(r'${\bf{(a)}\ Cube}\ Re=$'+qc["Re"], fontfamily='serif', loc='left', fontsize='medium')
    ax[0,1].set_title(r'${\bf{(b)}\ Sphere}\ Re=$'+qs["Re"], fontfamily='serif', loc='left', fontsize='medium')
    ax[1,0].set_title(r'${\bf{(c)}', fontfamily='serif', loc='left', fontsize='medium')
    ax[1,1].set_title(r'${\bf{(d)}', fontfamily='serif', loc='left', fontsize='medium')

    plt.subplots_adjust(wspace=0.05,hspace=0.1)
    savename+=formt
    plt.savefig(savename, bbox_inches='tight', dpi=qual,format=formt)
    print('Saving',savename);plt.close();plt.clf()


def plot_BF_cube_sphere_X(savename,qc,qs,lim):
    
    xmax = 1; xmin = -xmax
    ymax = 1; ymin = -ymax
    zz=64
    cmap = discrete_cmap(zz,'gray')
    
    fig,ax = plt.subplots(ncols=2, nrows=1)
    fig.set_size_inches(fig_width,fig_height/3.222)
    
    ax[0].axis('equal')
    ax[1].axis('equal')
      
    ax[0].set_ylim(ymin,ymax)
    ax[1].set_ylim(ymin,ymax)
        
    ax[0].set_xlim(xmin,xmax) 
    ax[1].set_xlim(xmin,xmax)
    
    ax[0].set_ylabel(r'$y$',labelpad=-0) 
    ax[0].set_xlabel(r'$z$',labelpad=-0)
    ax[1].set_xlabel(r'$z$',labelpad=-0)
    
    ax[1].set_yticklabels([])

    levels = np.linspace(lim[0], lim[1], zz)
    lvs = [0,1] # np.linspace(lim[0]*0.8, lim[1]*0.8,4)

    ##### PLOT A #####
    axp0 = ax[0].tricontourf(qc['x']['tr'], qc['x']['u'], cmap=cmap,levels=levels,extend='both')
    ax[0].tricontour(qc['x']['tr'], qc['x']['u'], levels=lvs,  colors='k', linewidths=0.2)
    ax[0].tricontour(qc['x']['tr'], qc['x']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes0 = inset_axes(ax[0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$u$', fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)

    axp1 = ax[1].tricontourf(qs['x']['tr'], qs['x']['u'], cmap=cmap,levels=levels,extend='both')
    ax[1].tricontour(qs['x']['tr'], qs['x']['u'], levels=lvs,  colors='k', linewidths=0.2)
    ax[1].tricontour(qs['x']['tr'], qs['x']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes1 = inset_axes(ax[1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$u$', fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
    
    ax[0].grid(True, c='gray', ls=':', lw=0.6)
    ax[1].grid(True, c='gray', ls=':', lw=0.6)
    
    ##### Patches #####
    pc1 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    ps1 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    
    ax[0].add_artist(pc1)
    ax[1].add_artist(ps1)
    
    # Move patches to the front
    pc1.set_zorder(2)
    ps1.set_zorder(2)
    
    ax[0].set_title(r'${\bf{(a)}\ Cube}\ Re=$'+qc["Re"], fontfamily='serif', loc='left', fontsize='medium')
    ax[1].set_title(r'${\bf{(b)}\ Sphere}\ Re=$'+qs["Re"], fontfamily='serif', loc='left', fontsize='medium')

    plt.subplots_adjust(wspace=0.05)
    savename+=formt
    plt.savefig(savename, bbox_inches='tight', dpi=qual,format=formt)
    print('Saving',savename);plt.close();plt.clf()

def pltBF_CS3(savename,qc,qs,fc,fs,lim):
    
    xmin = -2.; xmax = abs(xmin)*3
    ymin = -1.; ymax = 1.
    zz=32+1
    cmap = discrete_cmap(zz,'gray')
    
    fig,ax = plt.subplots(ncols=2, nrows=3)
    fig.set_size_inches(fig_width,fig_height)
    
    ax[0,0].axis('equal')
    ax[0,1].axis('equal')
    ax[1,0].axis('equal')
    ax[1,1].axis('equal')
    ax[2,0].axis('equal')
    ax[2,1].axis('equal')
      
    ax[0,0].set_ylim(ymin,ymax)
    ax[0,1].set_ylim(ymin,ymax)
    ax[1,0].set_ylim(ymin,ymax)
    ax[1,1].set_ylim(ymin,ymax)
    ax[2,0].set_ylim(ymin,ymax)
    ax[2,1].set_ylim(ymin,ymax)
    
    ax[0,0].set_xlim(xmin,xmax) 
    ax[0,1].set_xlim(xmin,xmax)
    ax[1,0].set_xlim(xmin,xmax) 
    ax[1,1].set_xlim(xmin,xmax)
    ax[2,0].set_xlim(xmin,xmax)
    ax[2,1].set_xlim(xmin,xmax)
    
    ax[0,0].set_ylabel(r'$y$',labelpad=-0) 
    ax[1,0].set_ylabel(r'$z$',labelpad=-0)
    ax[2,0].set_ylabel(r'$y$',labelpad=-0)
    ax[0,0].set_xlabel(r'$x$',labelpad=-0)
    ax[0,1].set_xlabel(r'$x$',labelpad=-0)
    ax[1,0].set_xlabel(r'$x$',labelpad=-0)
    ax[1,1].set_xlabel(r'$x$',labelpad=-0)
    ax[2,0].set_xlabel(r'$z$',labelpad=-0)
    ax[2,1].set_xlabel(r'$z$',labelpad=-0)
    
    ax[0,0].set_xticklabels([])
    ax[0,1].set_xticklabels([])
    
    ax[1,0].set_xticklabels([])
    ax[1,1].set_xticklabels([])

    ax[0,1].set_yticklabels([])
    ax[1,1].set_yticklabels([])
    ax[2,1].set_yticklabels([])

    levels = np.linspace(lim[0], lim[1], zz)
    lvs = [0,1] # np.linspace(lim[0]*0.8, lim[1]*0.8,4)

    ##### PLOT A #####
    axp0 = ax[0,0].tricontourf(qc['y']['tr'], qc['y'][fc[0]], cmap=cmap,levels=levels,extend='both')
    ax[0,0].tricontour(qc['y']['tr'], qc['y'][fc[0]],levels=[0],colors='k', linewidths=0.2, linestyles='dotted')
    ax[0,0].tricontour(qc['y']['tr'], qc['y'][fc[0]],levels=[1],colors='k', linewidths=0.2, linestyles='dashed')
    
    cbaxes0 = inset_axes(ax[0,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)

    axp1 = ax[1,0].tricontourf(qc['z']['tr'], qc['z'][fc[0]], cmap=cmap,levels=levels,extend='both')
    ax[1,0].tricontour(qc['z']['tr'], qc['z'][fc[0]], levels=lvs,  colors='k', linewidths=0.2)
    ax[1,0].tricontour(qc['z']['tr'], qc['z'][fc[0]], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes1 = inset_axes(ax[1,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
    
    ##### PLOT B #####
    axp2 = ax[0,1].tricontourf(qs['y']['tr'], qs['y'][fs[0]], cmap=cmap,levels=levels,extend='both')
    ax[0,1].tricontour(qs['y']['tr'], qs['y'][fs[0]], levels=lvs,  colors='k', linewidths=0.2)
    ax[0,1].tricontour(qs['y']['tr'], qs['y'][fs[0]], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes2 = inset_axes(ax[0,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp2,orientation='horizontal',cax=cbaxes2)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
    
    axp3 = ax[1,1].tricontourf(qs['z']['tr'], qs['z'][fs[0]], cmap=cmap,levels=levels,extend='both')
    ax[1,1].tricontour(qs['z']['tr'], qs['z'][fs[0]], levels=lvs,  colors='k', linewidths=0.2)
    ax[1,1].tricontour(qs['z']['tr'], qs['z'][fs[0]], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes3 = inset_axes(ax[1,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp3,orientation='horizontal',cax=cbaxes3)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-8)
    cbar.set_ticklabels(ticks)
    
    ##### PLOT C #####
    axp0 = ax[2,0].tricontourf(qc['x']['tr'], qc['x'][fc[0]], cmap=cmap,levels=levels,extend='both')
    ax[2,0].tricontour(qc['x']['tr'], qc['x'][fc[0]], levels=lvs,  colors='k', linewidths=0.2)
    ax[2,0].tricontour(qc['x']['tr'], qc['x'][fc[0]], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes0 = inset_axes(ax[2,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)

    axp1 = ax[2,1].tricontourf(qs['x']['tr'], qs['x'][fs[0]], cmap=cmap,levels=levels,extend='both')
    ax[2,1].tricontour(qs['x']['tr'], qs['x'][fs[0]], levels=lvs,  colors='k', linewidths=0.2)
    ax[2,1].tricontour(qs['x']['tr'], qs['x'][fs[0]], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes1 = inset_axes(ax[2,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
        
    ax[0,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[0,1].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,1].grid(True, c='gray', ls=':', lw=0.6)
    ax[2,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[2,1].grid(True, c='gray', ls=':', lw=0.6)
      
    ##### Patches #####
    pc1 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    pc2 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    pc3 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)

    ps1 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    ps2 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    ps3 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    
    ax[0,0].add_artist(pc1)
    ax[0,1].add_artist(ps1)
    ax[1,0].add_artist(pc2)
    ax[1,1].add_artist(ps2)
    ax[2,0].add_artist(pc3)
    ax[2,1].add_artist(ps3)
    
    # Move patches to the front
    pc1.set_zorder(2)
    ps1.set_zorder(2)
    pc2.set_zorder(2)
    ps2.set_zorder(2)
    pc3.set_zorder(2)
    ps3.set_zorder(2)
    
    ax[0,0].set_title(r'${\bf (a)}\ z=0 \quad Cube\ - Re=$'+qc["Re"], fontfamily='serif', loc='left', fontsize='medium')
    ax[0,1].set_title(r'${\bf (b)}\ z=0 \quad Sphere\ - Re=$'+qs["Re"], fontfamily='serif', loc='left', fontsize='medium')

    ax[1,0].set_title(r'${\bf (c)}\ y=0$', fontfamily='serif', loc='left', fontsize='medium')
    ax[1,1].set_title(r'${\bf (d)}\ y=0$', fontfamily='serif', loc='left', fontsize='medium')

    ax[2,0].set_title(r'${\bf (e)}\ x=0$', fontfamily='serif', loc='left', fontsize='medium')
    ax[2,1].set_title(r'${\bf (f)}\ x=0$', fontfamily='serif', loc='left', fontsize='medium')

    plt.subplots_adjust(wspace=0.05,hspace=0.22)
    savename+=formt
    plt.savefig(savename, bbox_inches='tight', dpi=qual,format=formt)
    print('Saving',savename);plt.close();plt.clf()


def plot_MODE_cube_sphere(savename,qc,qs,lim):
    
    xmin = -2.; xmax = 6. ; ymin = -2; ymax = 2
    xmin = -4.; xmax = 16. ; ymin = -3; ymax = 3
    zz=64+1
    cmap = discrete_cmap(zz,'seismic')
    
    fig,ax = plt.subplots(ncols=2, nrows=2)
    fig.set_size_inches(fig_width,fig_height/2.222)
    
    ax[0,0].axis('equal'); ax[0,0].set_xticklabels([])
    ax[0,1].axis('equal'); ax[0,1].set_yticklabels([])
    ax[1,0].axis('equal'); ax[0,1].set_xticklabels([])
    ax[1,1].axis('equal'); ax[1,1].set_yticklabels([])
    
    ax[0,0].set_ylim(ymin,ymax); ax[0,0].set_xlim(xmin,xmax)
    ax[0,1].set_ylim(ymin,ymax); ax[0,1].set_xlim(xmin,xmax)
    ax[1,0].set_ylim(ymin,ymax); ax[1,0].set_xlim(xmin,xmax)
    ax[1,1].set_ylim(ymin,ymax); ax[1,1].set_xlim(xmin,xmax)

    ax[0,0].set_ylabel(r'$y$',labelpad=-0) 
    ax[1,0].set_ylabel(r'$z$',labelpad=-0)
    ax[1,0].set_xlabel(r'$x$',labelpad=-0)
    ax[1,1].set_xlabel(r'$x$',labelpad=-0)
    
    levels = np.linspace(lim[0], lim[1], zz)
    lvs = [0,1] # np.linspace(lim[0]*0.8, lim[1]*0.8,4)

    ##### PLOT A #####
    axp0 = ax[0,0].tricontourf(qc['y']['tr'], qc['y']['w'], cmap=cmap,levels=levels,extend='both')
    #ax[0,0].tricontour(qc['y']['tr'], qc['y']['w'], levels=lvs,  colors='k', linewidths=0.2)
    #ax[0,0].tricontour(qc['y']['tr'], qc['y']['w'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes0 = inset_axes(ax[0,0], width="3%", height="26%", loc=2) # 1 upper right
    cbar = plt.colorbar(axp0,orientation='vertical',cax=cbaxes0)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$w$', fontsize=6,labelpad=-12,rotation=1)
    cbar.set_ticklabels(ticks)

    axp1 = ax[1,0].tricontourf(qc['z']['tr'], qc['z']['w'], cmap=cmap,levels=levels,extend='both')
    #ax[1,0].tricontour(qc['z']['tr'], qc['z']['w'], levels=lvs,  colors='k', linewidths=0.2)
    #ax[1,0].tricontour(qc['z']['tr'], qc['z']['w'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes1 = inset_axes(ax[1,0], width="3%", height="26%", loc=2)
    cbar = plt.colorbar(axp1,orientation='vertical',cax=cbaxes1)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$w$', fontsize=6,labelpad=-12,rotation=1)
    cbar.set_ticklabels(ticks)
    
    ##### PLOT B #####
    axp2 = ax[0,1].tricontourf(qs['y']['tr'], qs['y']['w'], cmap=cmap,levels=levels,extend='both')
    #ax[0,1].tricontour(qs['y']['tr'], qs['y']['u'], levels=lvs,  colors='k', linewidths=0.2)
    #ax[0,1].tricontour(qs['y']['tr'], qs['y']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes2 = inset_axes(ax[0,1], width="3%", height="26%", loc=2) 
    cbar = plt.colorbar(axp2,orientation='vertical',cax=cbaxes2)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$w$', fontsize=6,labelpad=-12,rotation=1)
    cbar.set_ticklabels(ticks)
    
    axp3 = ax[1,1].tricontourf(qs['z']['tr'], qs['z']['w'], cmap=cmap,levels=levels,extend='both')
    #ax[1,1].tricontour(qs['z']['tr'], qs['z']['u'], levels=lvs,  colors='k', linewidths=0.2)
    #ax[1,1].tricontour(qs['z']['tr'], qs['z']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes3 = inset_axes(ax[1,1], width="3%", height="26%", loc=2)
    cbar = plt.colorbar(axp3,orientation='vertical',cax=cbaxes3)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(r'$w$', fontsize=6,labelpad=-12,rotation=1)
    cbar.set_ticklabels(ticks)
    
    ax[0,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[0,1].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,1].grid(True, c='gray', ls=':', lw=0.6)
    
    ##### Patches #####
    pc1 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    pc2 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    ps1 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    ps2 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    
    ax[0,0].add_artist(pc1); ax[0,1].add_artist(ps1); pc1.set_zorder(2)
    ax[1,0].add_artist(pc2); ax[1,1].add_artist(ps2); ps1.set_zorder(2)
        
    ax[0,0].set_title(r'${\bf{(a)}\ Cube}\ Re=$'+qc["Re"], fontfamily='serif', loc='left', fontsize='medium')
    ax[0,1].set_title(r'${\bf{(b)}\ Sphere}\ Re=$'+qs["Re"], fontfamily='serif', loc='left', fontsize='medium')
    #ax[1,0].set_title(r'${\bf{(c)}\ Cube}\ Re=$', fontfamily='serif', loc='left', fontsize='medium')
    #ax[1,1].set_title(r'${\bf{(d)}\ Sphere}\ Re=$', fontfamily='serif', loc='left', fontsize='medium')

    plt.subplots_adjust(wspace=0.05,hspace=0.1)
    savename+=formt
    plt.savefig(savename, bbox_inches='tight', dpi=qual,format=formt)
    print('Saving',savename);plt.close();plt.clf()


def pltMD_CS3(savename,qc,qs,fc,fs,lim):
    
    xmin = -2.; xmax = abs(xmin)*3
    ymin = -1.; ymax = 1.
    zz=32+1 # zero is centered in the middle
    cmap = discrete_cmap(zz,'seismic')
    
    fig,ax = plt.subplots(ncols=2, nrows=3)
    fig.set_size_inches(fig_width,fig_height)
    
    ax[0,0].axis('equal')
    ax[0,1].axis('equal')
    ax[1,0].axis('equal')
    ax[1,1].axis('equal')
    ax[2,0].axis('equal')
    ax[2,1].axis('equal')
      
    ax[0,0].set_ylim(ymin,ymax)
    ax[0,1].set_ylim(ymin,ymax)
    ax[1,0].set_ylim(ymin,ymax)
    ax[1,1].set_ylim(ymin,ymax)
    ax[2,0].set_ylim(ymin,ymax)
    ax[2,1].set_ylim(ymin,ymax)
    
    ax[0,0].set_xlim(xmin,xmax) 
    ax[0,1].set_xlim(xmin,xmax)
    ax[1,0].set_xlim(xmin,xmax) 
    ax[1,1].set_xlim(xmin,xmax)
    ax[2,0].set_xlim(xmin,xmax)
    ax[2,1].set_xlim(xmin,xmax)
    
    ax[0,0].set_ylabel(r'$y$',labelpad=-0) 
    ax[1,0].set_ylabel(r'$z$',labelpad=-0)
    ax[2,0].set_ylabel(r'$y$',labelpad=-0)
    ax[0,0].set_xlabel(r'$x$',labelpad=-0)
    ax[0,1].set_xlabel(r'$x$',labelpad=-0)
    ax[1,0].set_xlabel(r'$x$',labelpad=-0)
    ax[1,1].set_xlabel(r'$x$',labelpad=-0)
    ax[2,0].set_xlabel(r'$z$',labelpad=-0)
    ax[2,1].set_xlabel(r'$z$',labelpad=-0)
    
    ax[0,0].set_xticklabels([])
    ax[0,1].set_xticklabels([])
    
    ax[1,0].set_xticklabels([])
    ax[1,1].set_xticklabels([])

    ax[0,1].set_yticklabels([])
    ax[1,1].set_yticklabels([])
    ax[2,1].set_yticklabels([])

    ##### PLOT A #####
    xmin = lim[0]*qc['y'][fc[0]].min()
    xmax = lim[1]*qc['y'][fc[0]].max()
    print(xmin,xmax)
    levels = np.linspace(xmin,xmax,zz)
    axp0 = ax[0,0].tricontourf(qc['y']['tr'], qc['y'][fc[0]], cmap=cmap,levels=levels,extend='both')
    
    cbaxes0 = inset_axes(ax[0,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    ticks = np.around([xmin,xmax], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)

    xmin = lim[0]*qc['z'][fc[0]].min()
    xmax = lim[1]*qc['z'][fc[0]].max()
    print(xmin,xmax)
    levels = np.linspace(xmin,xmax,zz)
    axp1 = ax[1,0].tricontourf(qc['z']['tr'], qc['z'][fc[0]], cmap=cmap,levels=levels,extend='both')
    
    cbaxes1 = inset_axes(ax[1,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    ticks = np.around([xmin,xmax], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
    
    ##### PLOT B #####
    xmin = lim[0]*qs['y'][fc[0]].min()
    xmax = lim[1]*qs['y'][fc[0]].max()
    print(xmin,xmax)
    levels = np.linspace(xmin,xmax,zz)
    axp2 = ax[0,1].tricontourf(qs['y']['tr'], qs['y'][fs[0]], cmap=cmap,levels=levels,extend='both')
    
    cbaxes2 = inset_axes(ax[0,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp2,orientation='horizontal',cax=cbaxes2)
    ticks = np.around([xmin,xmax], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
    
    xmin = lim[0]*qs['z'][fc[0]].min()
    xmax = lim[1]*qs['z'][fc[0]].max()
    levels = np.linspace(xmin,xmax,zz)
    axp3 = ax[1,1].tricontourf(qs['z']['tr'], qs['z'][fs[0]], cmap=cmap,levels=levels,extend='both')
    
    cbaxes3 = inset_axes(ax[1,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp3,orientation='horizontal',cax=cbaxes3)
    ticks = np.around([xmin,xmax], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-8)
    cbar.set_ticklabels(ticks)
    
    ##### PLOT C #####
    xmin = lim[0]*qc['x'][fc[0]].min()
    xmax = lim[1]*qc['x'][fc[0]].max()
    print(xmin,xmax)
    levels = np.linspace(xmin,xmax,zz)
    axp0 = ax[2,0].tricontourf(qc['x']['tr'], qc['x'][fc[0]], cmap=cmap,levels=levels,extend='both')
    
    cbaxes0 = inset_axes(ax[2,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    ticks = np.around([xmin,xmax], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)

    xmin = lim[0]*qs['x'][fc[0]].min()
    xmax = lim[1]*qs['x'][fc[0]].max()
    print(xmin,xmax)
    levels = np.linspace(xmin,xmax,zz)
    axp1 = ax[2,1].tricontourf(qs['x']['tr'], qs['x'][fs[0]], cmap=cmap,levels=levels,extend='both')
    
    cbaxes1 = inset_axes(ax[2,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    ticks = np.around([xmin,xmax], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
        
    ax[0,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[0,1].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,1].grid(True, c='gray', ls=':', lw=0.6)
    ax[2,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[2,1].grid(True, c='gray', ls=':', lw=0.6)
      
    ##### Patches #####
    pc1 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    pc2 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    pc3 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)

    ps1 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    ps2 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    ps3 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    
    ax[0,0].add_artist(pc1)
    ax[0,1].add_artist(ps1)
    ax[1,0].add_artist(pc2)
    ax[1,1].add_artist(ps2)
    ax[2,0].add_artist(pc3)
    ax[2,1].add_artist(ps3)
    
    # Move patches to the front
    pc1.set_zorder(2)
    ps1.set_zorder(2)
    pc2.set_zorder(2)
    ps2.set_zorder(2)
    pc3.set_zorder(2)
    ps3.set_zorder(2)
    
    ax[0,0].set_title(r'${\bf (a)}\ z=0 \quad Cube\ - Re=$'+qc["Re"], fontfamily='serif', loc='left', fontsize='medium')
    ax[0,1].set_title(r'${\bf (b)}\ z=0 \quad Sphere\ - Re=$'+qs["Re"], fontfamily='serif', loc='left', fontsize='medium')

    ax[1,0].set_title(r'${\bf (c)}\ y=0$', fontfamily='serif', loc='left', fontsize='medium')
    ax[1,1].set_title(r'${\bf (d)}\ y=0$', fontfamily='serif', loc='left', fontsize='medium')

    ax[2,0].set_title(r'${\bf (e)}\ x=0$', fontfamily='serif', loc='left', fontsize='medium')
    ax[2,1].set_title(r'${\bf (f)}\ x=0$', fontfamily='serif', loc='left', fontsize='medium')

    plt.subplots_adjust(wspace=0.05,hspace=0.22)
    savename+=formt
    plt.savefig(savename, bbox_inches='tight', dpi=qual,format=formt)
    print('Saving',savename);plt.close();plt.clf()
    
def pltBF_CS3(savename,qc,qs,fc,fs,lm, flip = False):
    
    xmin = -2.; xmax = abs(xmin)*3
    ymin = -1.; ymax = 1.
    zz=32+1
    cmap = discrete_cmap(zz,'gray')
    
    fig,ax = plt.subplots(ncols=2, nrows=3)
    fig.set_size_inches(fig_width,fig_height)
    
    ax[0,0].axis('equal')
    ax[0,1].axis('equal')
    ax[1,0].axis('equal')
    ax[1,1].axis('equal')
    ax[2,0].axis('equal')
    ax[2,1].axis('equal')
      
    ax[0,0].set_ylim(ymin,ymax)
    ax[0,1].set_ylim(ymin,ymax)
    ax[1,0].set_ylim(ymin,ymax)
    ax[1,1].set_ylim(ymin,ymax)
    ax[2,0].set_ylim(ymin,ymax)
    ax[2,1].set_ylim(ymin,ymax)
    
    ax[0,0].set_xlim(xmin,xmax) 
    ax[0,1].set_xlim(xmin,xmax)
    ax[1,0].set_xlim(xmin,xmax) 
    ax[1,1].set_xlim(xmin,xmax)
    ax[2,0].set_xlim(xmin,xmax)
    ax[2,1].set_xlim(xmin,xmax)
    
    ax[0,0].set_ylabel(r'$y$',labelpad=-0) 
    ax[1,0].set_ylabel(r'$z$',labelpad=-0)
    ax[2,0].set_ylabel(r'$y$',labelpad=-0)
    ax[0,0].set_xlabel(r'$x$',labelpad=-0)
    ax[0,1].set_xlabel(r'$x$',labelpad=-0)
    ax[1,0].set_xlabel(r'$x$',labelpad=-0)
    ax[1,1].set_xlabel(r'$x$',labelpad=-0)
    ax[2,0].set_xlabel(r'$z$',labelpad=-0)
    ax[2,1].set_xlabel(r'$z$',labelpad=-0)
    
    ax[0,0].set_xticklabels([])
    ax[0,1].set_xticklabels([])
    
    ax[1,0].set_xticklabels([])
    ax[1,1].set_xticklabels([])

    ax[0,1].set_yticklabels([])
    ax[1,1].set_yticklabels([])
    ax[2,1].set_yticklabels([])

    levels = np.linspace(lm[0], lm[1], zz)
    lvs = [0,1] # np.linspace(lim[0]*0.8, lim[1]*0.8,4)

    ##### PLOT A #####
    
    if flip:
        axp0 = ax[0,0].tricontourf(qc['z']['tr'], qc['z'][fc[0]], cmap=cmap,levels=levels,extend='both')
        ax[0,0].tricontour(qc['z']['tr'], qc['z'][fc[0]],levels=[0.0],colors='w', linewidths=0.5, linestyles='solid')
        ax[0,0].tricontour(qc['z']['tr'], qc['z'][fc[0]],levels=[1.0],colors='k', linewidths=0.5, linestyles='solid')
    else:
        axp0 = ax[0,0].tricontourf(qc['y']['tr'], qc['y'][fc[0]], cmap=cmap,levels=levels,extend='both')
        ax[0,0].tricontour(qc['y']['tr'], qc['y'][fc[0]],levels=[0.0],colors='w', linewidths=0.5, linestyles='solid')
        ax[0,0].tricontour(qc['y']['tr'], qc['y'][fc[0]],levels=[1.0],colors='k', linewidths=0.5, linestyles='solid')
        
    cbaxes0 = inset_axes(ax[0,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)

    if flip:
        axp1 = ax[1,0].tricontourf(qc['y']['tr'], qc['y'][fc[0]], cmap=cmap,levels=levels,extend='both')
        ax[1,0].tricontour(qc['y']['tr'], qc['y'][fc[0]],levels=[0.0],colors='w', linewidths=0.5, linestyles='solid')
        ax[1,0].tricontour(qc['y']['tr'], qc['y'][fc[0]],levels=[1.0],colors='k', linewidths=0.5, linestyles='solid')
    else:
        axp1 = ax[1,0].tricontourf(qc['z']['tr'], qc['z'][fc[0]], cmap=cmap,levels=levels,extend='both')
        ax[1,0].tricontour(qc['z']['tr'], qc['z'][fc[0]],levels=[0.0],colors='w', linewidths=0.5, linestyles='solid')
        ax[1,0].tricontour(qc['z']['tr'], qc['z'][fc[0]],levels=[1.0],colors='k', linewidths=0.5, linestyles='solid')
        
    cbaxes1 = inset_axes(ax[1,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
    
    ##### PLOT B #####
    axp2 = ax[0,1].tricontourf(qs['y']['tr'], qs['y'][fs[0]], cmap=cmap,levels=levels,extend='both')
    ax[0,1].tricontour(qs['y']['tr'], qs['y'][fs[0]],levels=[0.0],colors='w', linewidths=0.5, linestyles='solid')
    ax[0,1].tricontour(qs['y']['tr'], qs['y'][fs[0]],levels=[1.0],colors='k', linewidths=0.5, linestyles='solid')
    
    cbaxes2 = inset_axes(ax[0,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp2,orientation='horizontal',cax=cbaxes2)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
    
    axp3 = ax[1,1].tricontourf(qs['z']['tr'], qs['z'][fs[0]], cmap=cmap,levels=levels,extend='both')
    ax[1,1].tricontour(qs['z']['tr'], qs['z'][fs[0]],levels=[0.0],colors='w', linewidths=0.5, linestyles='solid')
    ax[1,1].tricontour(qs['z']['tr'], qs['z'][fs[0]],levels=[1.0],colors='k', linewidths=0.5, linestyles='solid')
    
    cbaxes3 = inset_axes(ax[1,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp3,orientation='horizontal',cax=cbaxes3)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-8)
    cbar.set_ticklabels(ticks)
    
    ##### PLOT C #####
    axp0 = ax[2,0].tricontourf(qc['x']['tr'], qc['x'][fc[0]], cmap=cmap,levels=levels,extend='both')
    ax[2,0].tricontour(qc['x']['tr'], qc['x'][fc[0]],levels=[0.0],colors='w', linewidths=0.5, linestyles='solid')
    ax[2,0].tricontour(qc['x']['tr'], qc['x'][fc[0]],levels=[1.0],colors='k', linewidths=0.5, linestyles='solid')
    
    cbaxes0 = inset_axes(ax[2,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)

    axp1 = ax[2,1].tricontourf(qs['x']['tr'], qs['x'][fs[0]], cmap=cmap,levels=levels,extend='both')
    ax[2,1].tricontour(qs['x']['tr'], qs['x'][fs[0]],levels=[0.0],colors='w', linewidths=0.5, linestyles='solid')
    ax[2,1].tricontour(qs['x']['tr'], qs['x'][fs[0]],levels=[1.0],colors='k', linewidths=0.5, linestyles='solid')
    
    cbaxes1 = inset_axes(ax[2,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
        
    ax[0,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[0,1].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,1].grid(True, c='gray', ls=':', lw=0.6)
    ax[2,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[2,1].grid(True, c='gray', ls=':', lw=0.6)
      
    ##### Patches #####
    pc1 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    pc2 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    pc3 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)

    ps1 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    ps2 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    ps3 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    
    ax[0,0].add_artist(pc1)
    ax[0,1].add_artist(ps1)
    ax[1,0].add_artist(pc2)
    ax[1,1].add_artist(ps2)
    ax[2,0].add_artist(pc3)
    ax[2,1].add_artist(ps3)
    
    # Move patches to the front
    pc1.set_zorder(2)
    ps1.set_zorder(2)
    pc2.set_zorder(2)
    ps2.set_zorder(2)
    pc3.set_zorder(2)
    ps3.set_zorder(2)
    
                        #  {(\textit{a})\ Cube}
    ax[0,0].set_title(r'(\textit{a})\ $z=0\ Cube\ - Re='+str(qc["Re"])+'$', fontfamily='serif', loc='left', fontsize='medium')
    ax[0,1].set_title(r'(\textit{b})\ $z=0\ Sphere\ - Re='+str(qs["Re"])+'$', fontfamily='serif', loc='left', fontsize='medium')
    ax[1,0].set_title(r'(\textit{c})\ $y=0$', fontfamily='serif', loc='left', fontsize='medium')
    ax[1,1].set_title(r'(\textit{d})\ $y=0$', fontfamily='serif', loc='left', fontsize='medium')
    ax[2,0].set_title(r'(\textit{e})\ $x=0$', fontfamily='serif', loc='left', fontsize='medium')
    ax[2,1].set_title(r'(\textit{f})\ $x=0$', fontfamily='serif', loc='left', fontsize='medium')
    plt.subplots_adjust(wspace=0.05,hspace=0.22)
    savename+=formt
    plt.savefig(savename, dpi=qual,format=formt)#bbox_inches='tight',
    print('Saving',savename);plt.close();plt.clf()

def pltMD_CS3(savename,qc,qs,fc,fs,lm):
    
    xmin = -2.; xmax = abs(xmin)*3
    ymin = -1.; ymax = 1.
    zz=32+1 # zero is centered in the middle
    cmap = discrete_cmap(zz,'seismic')
    
    fig,ax = plt.subplots(ncols=2, nrows=3)
    fig.set_size_inches(fig_width,fig_height)
    
    ax[0,0].axis('equal')
    ax[0,1].axis('equal')
    ax[1,0].axis('equal')
    ax[1,1].axis('equal')
    ax[2,0].axis('equal')
    ax[2,1].axis('equal')
      
    ax[0,0].set_ylim(ymin,ymax)
    ax[0,1].set_ylim(ymin,ymax)
    ax[1,0].set_ylim(ymin,ymax)
    ax[1,1].set_ylim(ymin,ymax)
    ax[2,0].set_ylim(ymin,ymax)
    ax[2,1].set_ylim(ymin,ymax)
    
    ax[0,0].set_xlim(xmin,xmax) 
    ax[0,1].set_xlim(xmin,xmax)
    ax[1,0].set_xlim(xmin,xmax) 
    ax[1,1].set_xlim(xmin,xmax)
    ax[2,0].set_xlim(xmin,xmax)
    ax[2,1].set_xlim(xmin,xmax)
    
    ax[0,0].set_ylabel(r'$y$',labelpad=-0) 
    ax[1,0].set_ylabel(r'$z$',labelpad=-0)
    ax[2,0].set_ylabel(r'$y$',labelpad=-0)
    ax[0,0].set_xlabel(r'$x$',labelpad=-0)
    ax[0,1].set_xlabel(r'$x$',labelpad=-0)
    ax[1,0].set_xlabel(r'$x$',labelpad=-0)
    ax[1,1].set_xlabel(r'$x$',labelpad=-0)
    ax[2,0].set_xlabel(r'$z$',labelpad=-0)
    ax[2,1].set_xlabel(r'$z$',labelpad=-0)
    
    ax[0,0].set_xticklabels([])
    ax[0,1].set_xticklabels([])
    
    ax[1,0].set_xticklabels([])
    ax[1,1].set_xticklabels([])

    ax[0,1].set_yticklabels([])
    ax[1,1].set_yticklabels([])
    ax[2,1].set_yticklabels([])

    levels = np.linspace(lm[0], lm[1], zz)
    lvs = [0.00000] # np.linspace(lim[0]*0.8, lim[1]*0.8,4)

    ##### PLOT A #####
    axp0 = ax[0,0].tricontourf(qc['y']['tr'], qc['y'][fc[0]], cmap=cmap,levels=levels,extend='both')
    #ax[0,0].tricontour(qc['y']['tr'], qc['y']['u'],levels=lvs            , colors='k', linewidths=0.2)
    #ax[0,0].tricontour(qc['y']['tr'], qc['y']['u'],levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes0 = inset_axes(ax[0,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)

    axp1 = ax[1,0].tricontourf(qc['z']['tr'], qc['z'][fc[0]], cmap=cmap,levels=levels,extend='both')
    #ax[1,0].tricontour(qc['z']['tr'], qc['z']['u'], levels=lvs,  colors='k', linewidths=0.2)
    #ax[1,0].tricontour(qc['z']['tr'], qc['z']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes1 = inset_axes(ax[1,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
    
    ##### PLOT B #####
    axp2 = ax[0,1].tricontourf(qs['y']['tr'], qs['y'][fs[0]], cmap=cmap,levels=levels,extend='both')
    #ax[0,1].tricontour(qs['y']['tr'], qs['y']['u'], levels=lvs,  colors='k', linewidths=0.2)
    #ax[0,1].tricontour(qs['y']['tr'], qs['y']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes2 = inset_axes(ax[0,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp2,orientation='horizontal',cax=cbaxes2)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
    
    axp3 = ax[1,1].tricontourf(qs['z']['tr'], qs['z'][fs[0]], cmap=cmap,levels=levels,extend='both')
    #ax[1,1].tricontour(qs['z']['tr'], qs['z']['u'], levels=lvs,  colors='k', linewidths=0.2)
    #ax[1,1].tricontour(qs['z']['tr'], qs['z']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes3 = inset_axes(ax[1,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp3,orientation='horizontal',cax=cbaxes3)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-8)
    cbar.set_ticklabels(ticks)
    
    ##### PLOT C #####
    axp0 = ax[2,0].tricontourf(qc['x']['tr'], qc['x'][fc[0]], cmap=cmap,levels=levels,extend='both')
    #ax[2,0].tricontour(qc['x']['tr'], qc['x']['u'], levels=lvs,  colors='k', linewidths=0.2)
    #ax[2,0].tricontour(qc['x']['tr'], qc['x']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes0 = inset_axes(ax[2,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)

    axp1 = ax[2,1].tricontourf(qs['x']['tr'], qs['x'][fs[0]], cmap=cmap,levels=levels,extend='both')
    #ax[2,1].tricontour(qs['x']['tr'], qs['x']['u'], levels=lvs,  colors='k', linewidths=0.2)
    #ax[2,1].tricontour(qs['x']['tr'], qs['x']['u'], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    cbaxes1 = inset_axes(ax[2,1], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    ticks = np.around([lm[0], lm[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(fs[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)
        
    ax[0,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[0,1].grid(True, c='gray', ls=':', lw=0.6)
    ax[1,1].grid(True, c='gray', ls=':', lw=0.6)
    ax[2,0].grid(True, c='gray', ls=':', lw=0.6)
    ax[2,1].grid(True, c='gray', ls=':', lw=0.6)
      
    ##### Patches #####
    pc1 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    pc2 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    pc3 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)

    ps1 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    ps2 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    ps3 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    
    ax[0,0].add_artist(pc1)
    ax[0,1].add_artist(ps1)
    ax[1,0].add_artist(pc2)
    ax[1,1].add_artist(ps2)
    ax[2,0].add_artist(pc3)
    ax[2,1].add_artist(ps3)
    
    # Move patches to the front
    pc1.set_zorder(2)
    ps1.set_zorder(2)
    pc2.set_zorder(2)
    ps2.set_zorder(2)
    pc3.set_zorder(2)
    ps3.set_zorder(2)
    
    ax[0,0].set_title(r'${\bf (a)}\ z=0 \quad Cube\ - Re='+str(qc["Re"])+'$', fontfamily='serif', loc='left', fontsize='medium')
    ax[0,1].set_title(r'${\bf (b)}\ z=0 \quad Sphere\ - Re='+str(qs["Re"])+'$', fontfamily='serif', loc='left', fontsize='medium')

    ax[1,0].set_title(r'${\bf (c)}\ y=0$', fontfamily='serif', loc='left', fontsize='medium')
    ax[1,1].set_title(r'${\bf (d)}\ y=0$', fontfamily='serif', loc='left', fontsize='medium')

    ax[2,0].set_title(r'${\bf (e)}\ x=0$', fontfamily='serif', loc='left', fontsize='medium')
    ax[2,1].set_title(r'${\bf (f)}\ x=0$', fontfamily='serif', loc='left', fontsize='medium')

    plt.subplots_adjust(wspace=0.05,hspace=0.22)
    savename+=formt
    plt.savefig(savename, bbox_inches='tight', dpi=qual,format=formt)
    print('Saving',savename);plt.close();plt.clf()





# def rdnek3d(fname, pl_pos=0.0, rtol=0, atol=1e-3, flip=False):
#     print('Reading file:', fname)
#     qo = readnek(fname)["data"]

#     # Isolate mesh in each direction
#     # print('x min max = ', xo.min(), xo.max())
#     # print('y min max = ', yo.min(), yo.max())
#     # print('z min max = ', zo.min(), zo.max())

#     # Isolate velocity components

#     # if flip:
#     #     print('Flipping v and w')
#     #     xo, zo, yo = qo[:,:,0].ravel(), qo[:,:,1].ravel(), qo[:,:,2].ravel()
#     #     uo, wo, vo = qo[:,:,3].ravel(), qo[:,:,4].ravel(), qo[:,:,5].ravel()
#     # else: # true -> change v with w (u stays the same)
#     xo, yo, zo = qo[:,:,0].ravel(), qo[:,:,1].ravel(), qo[:,:,2].ravel()
#     uo, vo, wo = qo[:,:,3].ravel(), qo[:,:,4].ravel(), qo[:,:,5].ravel()

#     # Correcting mesh coordinates close to zero
#     # tol = 1e-3
#     # xo = qo[:,:,0].ravel()
#     # xo[np.isclose(xo, 0, atol=tol)] = 0

#     # yo = qo[:,:,1].ravel()
#     # yo[np.isclose(yo, 0, atol=tol)] = 0

#     # zo = qo[:,:,2].ravel()
#     # zo[np.isclose(zo, 0, atol=tol)] = 0
    
#     # Apply Gaussian smoothing
#     #smoothing_sigma = 1
#     # xo = gaussian_filter(xo, smoothing_sigma)
#     # yo = gaussian_filter(yo, smoothing_sigma)
#     # zo = gaussian_filter(zo, smoothing_sigma)

#     # uo = gaussian_filter(uo, smoothing_sigma)
#     # vo = gaussian_filter(vo, smoothing_sigma)
#     # wo = gaussian_filter(wo, smoothing_sigma)
          
#     # print('')
#     # print('  u mn mx = ', uo.min(), uo.max())
#     # print(' v mn mx = ', vo.min(), vo.max())
#     # print('w mn mx = ', wo.min(), wo.max())
    
#     # Compute velocity magnitude
#     # U = np.sqrt(uo**2 + vo**2 + wo**2)
#     # print('U min max = ', U.min(), U.max())

#     q = {"x":{},"y":{},"x":{},"Re":{}} # initialize dictionary
    
#     plane = 'x'; #print('Plane', plane, '==', pl_pos)
#     inx = np.isclose(xo, pl_pos, rtol=rtol, atol=atol)
#     y, z, u, v, w = yo[inx], zo[inx], uo[inx], vo[inx], wo[inx]
#     tr = Triangulation(y, z, triangles=Delaunay(np.vstack((y, z)).T).simplices)
#     #tr = Triangulation(y, z, triangles=Delaunay(np.vstack((y, z)).T, qhull_options="QJ Qbb").simplices)
#     q[plane] = { 'y': y, 'z': z, 'u': u, 'v': v, 'w': w, 'tr': tr }

#     plane = 'y'
#     iny = np.isclose(yo, pl_pos, rtol=rtol, atol=atol)
#     x, z, u, v, w = xo[iny], zo[iny], uo[iny], vo[iny], wo[iny]
#     tr = Triangulation(x, z, triangles=Delaunay(np.vstack((x, z)).T).simplices)
#     q[plane] = { 'x': x, 'z': z, 'u': u, 'v': v, 'w': w, 'tr': tr }

#     plane = 'z'
#     inz = np.isclose(zo, pl_pos, rtol=rtol, atol=atol)
#     x, y, u, v, w = xo[inz], yo[inz], uo[inz], vo[inz], wo[inz]
#     tr = Triangulation(x, y, triangles=Delaunay(np.vstack((x, y)).T).simplices)
#     q[plane] = { 'x': x, 'y': y, 'u': u, 'v': v, 'w': w, 'tr': tr }
#     return q

def rdnek3d(fname, pl_pos=0.0, rtol=0, atol=1e-3):
    print('Reading file:', fname)
    qo = readnek(fname)["data"]

    xo, yo, zo = qo[:,:,0].ravel(), qo[:,:,1].ravel(), qo[:,:,2].ravel()
    uo, vo, wo = qo[:,:,3].ravel(), qo[:,:,4].ravel(), qo[:,:,5].ravel()

    # # Correcting mesh coordinates close to zero
    # if correct_zero:
    #     tol = 1e-2
    #     xo[np.isclose(xo, 0.0, atol=tol)] = 0.0
    #     yo[np.isclose(yo, 0.0, atol=tol)] = 0.0
    #     zo[np.isclose(zo, 0.0, atol=tol)] = 0.0

    # # Apply Gaussian smoothing
    # if apply_smoothing:
    #     xo = gaussian_filter(xo, smoothing_sigma)
    #     yo = gaussian_filter(yo, smoothing_sigma)
    #     zo = gaussian_filter(zo, smoothing_sigma)
    #     uo = gaussian_filter(uo, smoothing_sigma)
    #     vo = gaussian_filter(vo, smoothing_sigma)
    #     wo = gaussian_filter(wo, smoothing_sigma)

    q = {"x":{}, "y":{}, "z":{}, "Re":{}}

    for plane, axes in [('x', ('y', 'z')), ('y', ('x', 'z')), ('z', ('x', 'y'))]:
        axis1, axis2 = axes
        in_mask = np.isclose(qo[:,:,['x', 'y', 'z'].index(plane)].ravel(), pl_pos, rtol=rtol, atol=atol)

        # Directly index the 3D array for the specific plane data. Replace 3, 4, 5 with the correct indices for 'u', 'v', 'w'.
        axis1_data, axis2_data = [qo[:,:,['x', 'y', 'z'].index(a)].ravel()[in_mask] for a in [axis1, axis2]]
        u, v, w = [qo[:,:,i].ravel()[in_mask] for i in [3, 4, 5]]

        try:
            tr = Triangulation(axis1_data, axis2_data, triangles=Delaunay(np.vstack((axis1_data, axis2_data)).T).simplices)
            #tr = Triangulation(axis1_data, axis2_data, triangles=Delaunay(np.vstack((axis1_data, axis2_data)).T, qhull_options='Qt').simplices)
            #tr = Triangulation(axis1_data, axis2_data, triangles=Delaunay(np.vstack((axis1_data, axis2_data)).T, qhull_options='QJ').simplices)
        except Exception as e:
            print(f'Failed to create triangulation for plane {plane}: {e}')
            tr = None

        q[plane] = {axis1: axis1_data, axis2: axis2_data, 'u': u, 'v': v, 'w': w, 'tr': tr}

        # print(f'For plane {plane}: \n'
        #     f' u min, max = {u.min()}, {u.max()}\n'
        #     f' v min, max = {v.min()}, {v.max()}\n'
        #     f' w min, max = {w.min()}, {w.max()}')

    return q
