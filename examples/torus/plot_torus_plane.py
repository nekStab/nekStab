#!/usr/bin/env python3

#!/usr/bin/env python3
from bffunc import *

#pip install --upgrade ipython ipykernel ipympl matplotlib scipy django scienceplots
#conda isntall ipython ipykernel ipympl matplotlib numpy scipy django scienceplots

#https://github.com/matplotlib/matplotlib/issues/26287
# AttributeError: 'NoneType' object has no attribute '_get_renderer'
# bbox_inches='tight' commented

from time import time

def read_data(base_path, file_name, **kwargs):
    full_path = f"{base_path}/{file_name}"
    
    start_time = time()
    data = rdnek3d(full_path, **kwargs)
    elapsed_time = time() - start_time
    
    print(f"Elapsed time for rdnek3d: {elapsed_time}s")
    data["Re"] = base_path.split('/')[-1]
    
    return data

def rdnek3d(fname, pl_pos=0.0, rtol=0.1, atol=1e-3):
    print('Reading file:', fname)
    qo = readnek(fname)["data"]
    xo, yo, zo = qo[:,:,0].ravel(), qo[:,:,1].ravel(), qo[:,:,2].ravel()
    uo, vo, wo = qo[:,:,3].ravel(), qo[:,:,4].ravel(), qo[:,:,5].ravel()
    q = {"x":{}, "y":{}, "z":{}, "Re":{}}

    for plane, axes in [('x', ('y', 'z'))]: # , ('y', ('x', 'z')), ('z', ('x', 'y'))
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

        print(f'For plane {plane}: \n'
            f' u min, max = {u.min()}, {u.max()}\n'
            f' v min, max = {v.min()}, {v.max()}\n'
            f' w min, max = {w.min()}, {w.max()}')

    return q

def pltBF_CS3(savename,qc,qs,fc,fs,lim):
    
    xmin = 2.5; xmax = 3.5
    ymin = -1.; ymax = 1.
    zz=256+1
    cmap = discrete_cmap(zz,'viridis')
    
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
    axp0 = ax[0,0].tricontourf(qc['x']['tr'], qc['x'][fc[0]], cmap=cmap,levels=levels,extend='both')
    ax[0,0].tricontour(qc['x']['tr'], qc['x'][fc[0]],levels=[0],colors='k', linewidths=0.2, linestyles='dotted')
    ax[0,0].tricontour(qc['x']['tr'], qc['x'][fc[0]],levels=[1],colors='k', linewidths=0.2, linestyles='dashed')
    
    cbaxes0 = inset_axes(ax[0,0], width="18%", height="7%", loc=2) 
    cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    ticks = np.around([lim[0], lim[1]], decimals=2)
    cbar.set_ticks(ticks)
    cbar.ax.tick_params(labelsize=6)
    # cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    cbar.set_ticklabels(ticks)

    # axp1 = ax[1,0].tricontourf(qc['x']['tr'], qc['x'][fc[0]], cmap=cmap,levels=levels,extend='both')
    # ax[1,0].tricontour(qc['x']['tr'], qc['x'][fc[0]], levels=lvs,  colors='k', linewidths=0.2)
    # ax[1,0].tricontour(qc['x']['tr'], qc['x'][fc[0]], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    # cbaxes1 = inset_axes(ax[1,0], width="18%", height="7%", loc=2) 
    # cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    # ticks = np.around([lim[0], lim[1]], decimals=2)
    # cbar.set_ticks(ticks)
    # cbar.ax.tick_params(labelsize=6)
    # cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    # cbar.set_ticklabels(ticks)
    
    # ##### PLOT B #####
    # axp2 = ax[0,1].tricontourf(qs['y']['tr'], qs['y'][fs[0]], cmap=cmap,levels=levels,extend='both')
    # ax[0,1].tricontour(qs['y']['tr'], qs['y'][fs[0]], levels=lvs,  colors='k', linewidths=0.2)
    # ax[0,1].tricontour(qs['y']['tr'], qs['y'][fs[0]], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    # cbaxes2 = inset_axes(ax[0,1], width="18%", height="7%", loc=2) 
    # cbar = plt.colorbar(axp2,orientation='horizontal',cax=cbaxes2)
    # ticks = np.around([lim[0], lim[1]], decimals=2)
    # cbar.set_ticks(ticks)
    # cbar.ax.tick_params(labelsize=6)
    # cbar.set_label(fs[1], fontsize=7,labelpad=-11)
    # cbar.set_ticklabels(ticks)
    
    # axp3 = ax[1,1].tricontourf(qs['z']['tr'], qs['z'][fs[0]], cmap=cmap,levels=levels,extend='both')
    # ax[1,1].tricontour(qs['z']['tr'], qs['z'][fs[0]], levels=lvs,  colors='k', linewidths=0.2)
    # ax[1,1].tricontour(qs['z']['tr'], qs['z'][fs[0]], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    # cbaxes3 = inset_axes(ax[1,1], width="18%", height="7%", loc=2) 
    # cbar = plt.colorbar(axp3,orientation='horizontal',cax=cbaxes3)
    # ticks = np.around([lim[0], lim[1]], decimals=2)
    # cbar.set_ticks(ticks)
    # cbar.ax.tick_params(labelsize=6)
    # cbar.set_label(fs[1], fontsize=7,labelpad=-8)
    # cbar.set_ticklabels(ticks)
    
    # ##### PLOT C #####
    # axp0 = ax[2,0].tricontourf(qc['x']['tr'], qc['x'][fc[0]], cmap=cmap,levels=levels,extend='both')
    # ax[2,0].tricontour(qc['x']['tr'], qc['x'][fc[0]], levels=lvs,  colors='k', linewidths=0.2)
    # ax[2,0].tricontour(qc['x']['tr'], qc['x'][fc[0]], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    # cbaxes0 = inset_axes(ax[2,0], width="18%", height="7%", loc=2) 
    # cbar = plt.colorbar(axp0,orientation='horizontal',cax=cbaxes0)
    # ticks = np.around([lim[0], lim[1]], decimals=2)
    # cbar.set_ticks(ticks)
    # cbar.ax.tick_params(labelsize=6)
    # cbar.set_label(fc[1], fontsize=7,labelpad=-11)
    # cbar.set_ticklabels(ticks)

    # axp1 = ax[2,1].tricontourf(qs['x']['tr'], qs['x'][fs[0]], cmap=cmap,levels=levels,extend='both')
    # ax[2,1].tricontour(qs['x']['tr'], qs['x'][fs[0]], levels=lvs,  colors='k', linewidths=0.2)
    # ax[2,1].tricontour(qs['x']['tr'], qs['x'][fs[0]], levels=-np.flipud(lvs), colors='k', linewidths=0.2)
    
    # cbaxes1 = inset_axes(ax[2,1], width="18%", height="7%", loc=2) 
    # cbar = plt.colorbar(axp1,orientation='horizontal',cax=cbaxes1)
    # ticks = np.around([lim[0], lim[1]], decimals=2)
    # cbar.set_ticks(ticks)
    # cbar.ax.tick_params(labelsize=6)
    # cbar.set_label(fs[1], fontsize=7,labelpad=-11)
    # cbar.set_ticklabels(ticks)
        
    # ax[0,0].grid(True, c='gray', ls=':', lw=0.6)
    # ax[1,0].grid(True, c='gray', ls=':', lw=0.6)
    # ax[0,1].grid(True, c='gray', ls=':', lw=0.6)
    # ax[1,1].grid(True, c='gray', ls=':', lw=0.6)
    # ax[2,0].grid(True, c='gray', ls=':', lw=0.6)
    # ax[2,1].grid(True, c='gray', ls=':', lw=0.6)
      
    # ##### Patches #####
    # pc1 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    # pc2 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)
    # pc3 = mpl.patches.Rectangle((-0.5,-0.5), 1,1, facecolor='gray', edgecolor='black', lw=0.6)

    # ps1 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    # ps2 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    # ps3 = mpl.patches.Circle((0,0), 0.5, facecolor='gray', edgecolor='black', lw=0.6)
    
    # ax[0,0].add_artist(pc1)
    # ax[0,1].add_artist(ps1)
    # ax[1,0].add_artist(pc2)
    # ax[1,1].add_artist(ps2)
    # ax[2,0].add_artist(pc3)
    # ax[2,1].add_artist(ps3)
    
    # # Move patches to the front
    # pc1.set_zorder(2)
    # ps1.set_zorder(2)
    # pc2.set_zorder(2)
    # ps2.set_zorder(2)
    # pc3.set_zorder(2)
    # ps3.set_zorder(2)
    
    ax[0,0].set_title(r'${\bf (a)}\ x=0 \quad - Re=$'+qc["Re"], fontfamily='serif', loc='left', fontsize='medium')
    ax[0,1].set_title(r'${\bf (b)}\ x=0 \quad - Re=$'+qs["Re"], fontfamily='serif', loc='left', fontsize='medium')

    ax[1,0].set_title(r'${\bf (c)}\ x=0$', fontfamily='serif', loc='left', fontsize='medium')
    ax[1,1].set_title(r'${\bf (d)}\ x=0$', fontfamily='serif', loc='left', fontsize='medium')

    ax[2,0].set_title(r'${\bf (e)}\ x=0$', fontfamily='serif', loc='left', fontsize='medium')
    ax[2,1].set_title(r'${\bf (f)}\ x=0$', fontfamily='serif', loc='left', fontsize='medium')

    plt.subplots_adjust(wspace=0.05,hspace=0.22)
    savename+=formt
    plt.savefig(savename, bbox_inches='tight', dpi=qual,format=formt)
    print('Saving',savename);plt.close();plt.clf()



q = read_data('.','BF_torus0.f00001', pl_pos=0.0, rtol=0.0, atol=1e-6)

q2 = read_data('.','dRetorus0.f00001', pl_pos=0.0, rtol=0.0, atol=1e-6)

pltBF_CS3('BF_torus.',q,q2,'u','v',[0,1])