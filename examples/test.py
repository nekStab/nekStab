#!/usr/bin/env python
import sys, os, stat, re, time
from unittest import skip
from shutil import copyfile
import unittest, argparse, shutil, fileinput
from subprocess import call, PIPE, STDOUT, Popen
from subprocess import check_call as ccall 

def rnek(cwd, par_file, ifmpi, log_suffix="", n_procs=1, step_limit=None, verbose=False):
    toc = time.perf_counter()
    try:
        nek5000 = os.path.join(cwd, "nek5000")
        logfile = os.path.join(cwd, 'logfile.{0}{1}'.format(n_procs, log_suffix))
        compfile = os.path.join(cwd, "compfile")
        session_name = os.path.join(cwd, "SESSION.NAME")
        ioinfo = os.path.join(cwd, "ioinfo")
        print("Compiling nek5000 at", cwd,'> compfile')
        with open(compfile, "w") as f1:
            command = ["./cmpile.sh", "all"]
            call(command, cwd=cwd, stdout=f1)
        if ifmpi:
            command = ["mpiexec", "-np", str(n_procs), nek5000]
        else:
            command = [nek5000]
        print("Running nek5000...")
        print('    Using command "{0}"'.format(" ".join(command)))
        print('    Using working directory "{0}"'.format(cwd))
        print('    Using file "{0}"'.format(par_file))
        #p_file = os.path.join(cwd,par_file+'.par')
        #bp_file = p_file+'.{0}{1}'.format(n_procs, log_suffix)
        #shutil.copyfile(p_file,bp_file)
        #print('    Backup to file "{0}"'.format(bp_file))
        print()
        print(' tail -f',logfile)
        try:
            with open(session_name, "w") as f:
                f.writelines(
                    [
                    "{0}\n".format(1),
                    "{0}\n".format(par_file),
                    "{0}\n".format(cwd + "/"),
                    ]
                )

            if step_limit:
                with open(ioinfo, "w") as f:
                    f.writelines(["-{0}".format(step_limit)])
    
            if verbose:
                with open(logfile, "w") as f:
                    proc = Popen(command, cwd=cwd, stderr=STDOUT, stdout=PIPE)
                    for line in proc.stdout:
                        sys.stdout.write(line)
                        f.write(line)
            else:
                with open(logfile, "w") as f:
                    call(command, cwd=cwd, stdout=f)

        except Exception as E:
            # TODO: Change to warnings.warn()
            print(" ###### Error: {0}".format(E))
            print()
            print()
            print()
            print()
        else:
            toc = time.perf_counter(); ttime=toc-tic
            print(f"Case finished in in {ttime:0.2f} seconds")
            print(f"                      {ttime/60:0.1f} minutes")
            print(f"                      {ttime/3600:0.2f} hours")
            print('OK')
            print('OK')
            print('OK')
    except:
        pass
###############################################################
if __name__ == "__main__":
    
    tic = time.perf_counter()

    print('nohup python3 -u full_test.py >>logfile 2>&1 &')
    
    root = os.getcwd()  # main folder loation
    print("Current working directory: {0}".format(os.getcwd()))

    nps = 6
    
    # DNS
    rnek(root+'/cylinder/dns','1cyl',True,n_procs=nps)

    # BASE FLOW
    rnek(root+'/cylinder/baseflow/sfd'       ,'1cyl',True,n_procs=nps)
    rnek(root+'/cylinder/baseflow/sfd_dyn'   ,'1cyl',True,n_procs=nps)
    rnek(root+'/cylinder/baseflow/sfd_oifs'  ,'1cyl',True,n_procs=nps)
    rnek(root+'/cylinder/baseflow/boostconv' ,'1cyl',True,n_procs=nps)
    rnek(root+'/cylinder/baseflow/newton'    ,'1cyl',True,n_procs=nps)
    rnek(root+'/cylinder/baseflow/newton_dyn','1cyl',True,n_procs=nps)
    rnek(root+'/cylinder/baseflow/newton_upo','1cyl',True,n_procs=nps)

    # STABILITY 
    rnek(root+'/cylinder/stability/direct'         ,'1cyl',True,n_procs=nps)
    rnek(root+'/cylinder/stability/adjoint'        ,'1cyl',True,n_procs=nps)
    rnek(root+'/cylinder/stability/direct_Floquet' ,'1cyl',True,n_procs=nps)
    rnek(root+'/cylinder/stability/adjoint_Floquet','1cyl',True,n_procs=nps)

    # BUDGET,WAVEMAKER,BASEFLOW SENSITIVITY
    #rnek(root+'/4postprocessing/sensitivity',cn,True,n_procs=nps)
    #rnek(root+'/4postprocessing/steady_force_sensitivity',cn,True,n_procs=nps)
    #rnek(root+'/4postprocessing/steady_force_sensitivity',cn,True,n_procs=nps)
   
    # BACKWARD FACING STEP
    rnek(root+'/back_fstep/baseflow'        ,'bfs',True,n_procs=nps)
    rnek(root+'/back_fstep/transient_growth','bfs',True,n_procs=nps)
    
    # THERMOSYPHON
    rnek(root+'/thersyphon/baseflow','tsyphon',True,n_procs=nps)
    
    # AXISYMETRIC JET
    rnek(root+'/tpjet/baseflow'     ,'tpjet',True,n_procs=nps)

    # 2 CYLINDERS
    rnek(root+'/flip_flop/baseflow' ,'2cyl',True,n_procs=nps)

    toc = time.perf_counter(); ttime=toc-tic
    print(f"Script finished in in {ttime:0.2f} seconds")
    print(f"                      {ttime/60:0.1f} minutes")
    print(f"                      {ttime/3600:0.2f} hours")
