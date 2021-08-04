#!/usr/bin/env python
import sys, os, stat, re, time
from unittest import skip
from shutil import copyfile
import unittest, argparse, shutil, fileinput
from subprocess import call, PIPE, STDOUT, Popen
from subprocess import check_call as ccall 

def rnek(cwd, par_file, ifmpi, log_suffix="", n_procs=1, step_limit=None, verbose=False):
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
        print("Could not successfully run nek5000! Caught error: {0}".format(E))
    else:
        print("Finished running nek5000!")

###############################################################
if __name__ == "__main__":
    
    tic = time.perf_counter()

    print('nohup python3 -u full_test.py >>logfile 2>&1 &')
    
    root = os.getcwd()  # main folder loation
    print("Current working directory: {0}".format(os.getcwd()))

    cn = "1cyl"  # case name
    nps = 6
    
    # DNS
    rnek(root+'/0dns/',cn,True,n_procs=nps)

    #SFD
    rnek(root+'/1_2baseflow/sfd',cn,True,n_procs=nps)

    #SFD OIFS 
    rnek(root+'/1_2baseflow/sfd_oifs',cn,True,n_procs=nps)
    
    #SFD OIFS + adaptive
    rnek(root+'/1_2baseflow/sfd_dyn',cn,True,n_procs=nps)

    # BOOSTCONV
    rnek(root+'/1_2baseflow/boostconv',cn,True,n_procs=nps)

    # NEWTON dynamic
    rnek(root+'/1_2baseflow/newton_dyn',cn,True,n_procs=nps)
    
    # NEWTON
    rnek(root+'/1_2baseflow/newton',cn,True,n_procs=nps)

    # DIRECT 
    rnek(root+'/3stability/direct',cn,True,n_procs=nps)

    # ADJOINT
    rnek(root+'/3stability/adjoint',cn,True,n_procs=nps)

    # WAVEMAKER
    rnek(root+'/4postprocessing/wavemaker',cn,True,n_procs=nps)

    # BUDGET
    rnek(root+'/4postprocessing/energy_budget',cn,True,n_procs=nps)
        
    # BF SENSITIVITY
    rnek(root+'/4postprocessing/baseflow_sensitivity',cn,True,n_procs=nps)

    # FORCING SENSITIVITY
    rnek(root+'/4postprocessing/steady_force_sensitivity',cn,True,n_procs=nps)

    toc = time.perf_counter(); ttime=toc-tic
    print(f"Script finished in in {ttime:0.2f} seconds")
    print(f"                      {ttime/60:0.1f} minutes")
    print(f"                      {ttime/3600:0.2f} hours")