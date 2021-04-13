#!/usr/bin/env python
import sys, os, stat, re, time
from unittest import skip
from shutil import copyfile
import unittest, argparse, shutil, fileinput
from subprocess import call, PIPE, STDOUT, Popen
from subprocess import check_call as ccall 
from numpy import exp, log10

CEND      = '\33[0m'
CBOLD     = '\33[1m'
CITALIC   = '\33[3m'
CURL      = '\33[4m'
CBLINK    = '\33[5m'
CBLINK2   = '\33[6m'
CSELECTED = '\33[7m'

CBLACK  = '\33[30m'
CRED    = '\33[31m'
CGREEN  = '\33[32m'
CYELLOW = '\33[33m'
CBLUE   = '\33[34m'
CVIOLET = '\33[35m'
CBEIGE  = '\33[36m'
CWHITE  = '\33[37m'

def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.exists(d):
            try:
                shutil.rmtree(d)
            except Exception as e:
                print(e)
                os.unlink(d)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)
    #shutil.rmtree(src)

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
    p_file = os.path.join(cwd,par_file+'.par')
    bp_file = p_file+'.{0}{1}'.format(n_procs, log_suffix)
    shutil.copyfile(p_file,bp_file)
    print('    Backup to file "{0}"'.format(bp_file))
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


def cSZ(infile, outfile,params):
    with open(infile, "r") as f:
        lines = f.readlines()
    # Substitute all the variables
    for key, value in params.items():
        print("Modf.",key, "to",value,'in',infile)
        if value:
            lines = [
                re.sub(
                    r"(.*\bparameter\b.*\b{0} *= *)\S+?( *[),])".format(key),
                    r"\g<1>{0}\g<2>".format(value),
                    l,
                    flags=re.I,
                )
                for l in lines
            ]
    with open(outfile, "w") as f:
        f.writelines(lines)

def c_pf(infile, outfile, opts):
    import configparser
    parfile = configparser.ConfigParser()
    parfile.read(infile)
    for section, name_vals in opts.items():
        for name, val in name_vals.items():
            print("In",section,':',name,'set to',val,'in',infile)
            parfile.set(section, name, val)
    with open(outfile, "w") as f:
        parfile.write(f)

def ccopy(base,path):
    try:
        shutil.copytree(base, path)
        #copytree(base, path)
        print("Copying", base, "to", path)
    except OSError:
        pass

def rec(x):
        return 3.83207915e+03 -6.96661369e+02*x  +5.99070559e+01*x**2 -3.00364213e+00*x**3  +9.47417747e-02*x**4 -1.93536665e-03*x**5 +2.55607403e-05*x**6 -2.10608078e-07*x**7 + 9.84008556e-10*x**8 -1.99034184e-12*x**9
def Stt(x):
        return 2.65677490e+00 -4.11598506e-01*x +  3.52670003e-02*x**2 -1.83464151e-03*x**3 + 6.05343182e-05*x**4 -1.28872800e-06*x**5 +1.75989973e-08*x**6 -1.48597380e-10*x**7 + 7.05235396e-13*x**8 -1.43724492e-15*x**9 

###############################################################
if __name__ == "__main__":
    tic = time.perf_counter()

    print('to run this script:')
    print('nohup python3 -u autobot_naca.py >>logfile 2>&1 &')
    root = os.getcwd()  # main folder loation
    print("Current working directory: {0}".format(os.getcwd()))
    base = "base"  # reference case - base case to copy
    cn = "naca4412"  # case name

    aoa = ['8']#,'10','12','14','16','18','20','24','28','32','40','50','60']
    nps = 6

    for a in range(len(aoa)):
        estimated_rec = round(rec(float(aoa[a])))
        print('AoA=',aoa[a],'Re_c=',estimated_rec)
        rey = [round(0.65*estimated_rec),round(0.85*estimated_rec)]
        print('Re_sub=',rey[0])
        print('Re_sup=',rey[1])

        for r in range(len(rey)):
            ctic = time.perf_counter()

            path = 'aoa_'+str(aoa[a]) + "_Re_" + str(rey[r])
            ccopy(base,path)
            
            folder = root +'/' + path + '/'
            pf = folder + cn + ".par"
            SZ = folder + "SIZE"
            # {section: {optname: value, ...}, ...}
            print(pf)

            # PARAMS
            Re = rey[r]
            Pe = round(0.71*Re)
            St = round(abs(Stt(float(aoa[a]))),4)
            Tau = round((1/St),2)
            print('Re=',Re,' Pe=',Pe,' St=',St,' Tau=',Tau)

            shutil.copyfile('mshs/2d_a'+str(aoa[a])+'/naca4412.re2',folder+'/naca4412.re2')
            shutil.copyfile('mshs/2d_a'+str(aoa[a])+'/naca4412.ma2',folder+'/naca4412.ma2')

            # DNS
            tole = '1.E-7'; tolep = '1e-5'
            if r == 0: # subcritical -- start from zero
                c_pf(pf,pf,{'GENERAL':{'startFrom':'0'}})
            else: # supercritical -- start from previous case
                prev_path = 'aoA_'+str(aoa[a]) + "_Re_" + str(rey[r-1])
                prev_folder = root +'/' + prev_path + '/'
                shutil.copyfile(prev_folder+'/naca44120.f00030',folder+'/naca44120.f00001')
                c_pf(pf,pf,{'GENERAL':{'startFrom':'naca44120.f00001'}})

            c_pf(pf,pf,{'GENERAL':{'endTime':str(round(30*Tau,2))}})
            c_pf(pf,pf,{'GENERAL':{'writeInterval':str(round(Tau,2))}})
            c_pf(pf,pf,{'GENERAL':{'userParam01':'0'}})
            c_pf(pf,pf,{'GENERAL':{'userParam03':'0'}})
            c_pf(pf,pf,{'GENERAL':{'variableDt':'yes'}})
            c_pf(pf,pf,{'GENERAL':{'targetCFL':'1.5'}})
            c_pf(pf,pf,{'GENERAL':{'extrapolation':'OIFS'}})
            c_pf(pf,pf,{'VELOCITY':{'viscosity':'-'+str(Re)}})
            c_pf(pf,pf,{'VELOCITY':{'residualtol':tole}})
            c_pf(pf,pf,{'VELOCITY':{'residualproj':'no'}})
            c_pf(pf,pf,{'PRESSURE':{'residualtol':tolep}})
            c_pf(pf,pf,{'PRESSURE':{'residualproj':'no'}})
            c_pf(pf,pf,{'TEMPERATURE':{'conductivity':'-'+str(Pe)}})
            c_pf(pf,pf,{'TEMPERATURE':{'residualtol':tole}})
            c_pf(pf,pf,{'TEMPERATURE':{'residualproj':'no'}})
            rnek(folder,cn,True,log_suffix="_0dns",n_procs=nps)
            ccall(('visnek '+cn).split(), cwd=folder)
            
            # BASE FLOW
            tole = '1.E-10'
            # BASE FLOW # sfd
            c_pf(pf,pf,{'GENERAL':{'endTime':'10000'}})
            c_pf(pf,pf,{'GENERAL':{'userParam01':'1'}})
            c_pf(pf,pf,{'GENERAL':{'userParam03':'1'}})
            c_pf(pf,pf,{'GENERAL':{'userParam04':str(Tau)}})
            c_pf(pf,pf,{'GENERAL':{'userParam05':'0.1'}})
            c_pf(pf,pf,{'GENERAL':{'variableDt':'yes'}})
            c_pf(pf,pf,{'GENERAL':{'targetCFL':'1.5'}})
            c_pf(pf,pf,{'GENERAL':{'extrapolation':'OIFS'}})
            c_pf(pf,pf,{'GENERAL':{'writeInterval':'9999'}})
            c_pf(pf,pf,{'VELOCITY':{'viscosity':'-'+str(Re)}})
            c_pf(pf,pf,{'VELOCITY':{'residualtol':tole}})
            c_pf(pf,pf,{'VELOCITY':{'residualproj':'yes'}})
            c_pf(pf,pf,{'PRESSURE':{'residualtol':tole}})
            c_pf(pf,pf,{'PRESSURE':{'residualproj':'yes'}})
            c_pf(pf,pf,{'TEMPERATURE':{'conductivity':'-'+str(Pe)}})
            c_pf(pf,pf,{'TEMPERATURE':{'residualtol':tole}})
            c_pf(pf,pf,{'TEMPERATURE':{'residualproj':'yes'}})
            rnek(folder,cn,True,log_suffix="_1b",n_procs=nps)
            ccall(('visnek BF_'+cn).split(), cwd=folder)

            # DIRECT 
            tole = tole ; tolep = '1.0E-7'
            c_pf(pf,pf,{'GENERAL':{'endTime':str(Tau/8)}}) # T/8
            c_pf(pf,pf,{'GENERAL':{'userParam01':'3'}})
            c_pf(pf,pf,{'GENERAL':{'userParam03':'0'}})
            c_pf(pf,pf,{'GENERAL':{'variableDt':'no'}})
            c_pf(pf,pf,{'GENERAL':{'targetCFL':'0.5'}})
            c_pf(pf,pf,{'GENERAL':{'extrapolation':'standard'}})
            c_pf(pf,pf,{'VELOCITY':{'residualtol':tole}})
            c_pf(pf,pf,{'VELOCITY':{'residualproj':'yes'}})
            c_pf(pf,pf,{'PRESSURE':{'residualtol':tolep}})
            c_pf(pf,pf,{'PRESSURE':{'residualproj':'yes'}})
            c_pf(pf,pf,{'TEMPERATURE':{'residualtol':tole}})
            c_pf(pf,pf,{'TEMPERATURE':{'residualproj':'yes'}})
            rnek(folder,cn,True,log_suffix="_2d",n_procs=nps)
            ccall(('visnek dRe'+cn).split(), cwd=folder)
            ccall(('visnek dIm'+cn).split(), cwd=folder)

            # ADJOINT
            c_pf(pf,pf,{'GENERAL':{'userParam01':'3.2'}})
            rnek(folder,cn,True,log_suffix="_3a",n_procs=nps)
            ccall(('visnek aRe'+cn).split(), cwd=folder)
            ccall(('visnek aIm'+cn).split(), cwd=folder)

            # WAVEMAKER
            c_pf(pf,pf,{'GENERAL':{'userParam01':'4.1'}})
            rnek(folder,cn,True,log_suffix="_4w",n_procs=nps)
            ccall(('visnek wm_'+cn).split(), cwd=folder)

            # BF SENSITIVITY
            c_pf(pf,pf,{'GENERAL':{'userParam01':'4.2'}})
            rnek(folder,cn,True,log_suffix="_5s",n_procs=nps)
            ccall(('visnek tr_'+cn).split(), cwd=folder)
            ccall(('visnek ti_'+cn).split(), cwd=folder)
            ccall(('visnek pr_'+cn).split(), cwd=folder)
            ccall(('visnek pi_'+cn).split(), cwd=folder)
            ccall(('visnek sr_'+cn).split(), cwd=folder)
            ccall(('visnek si_'+cn).split(), cwd=folder)
            
            # FORCING SENSITIVITY

            c_pf(pf,pf,{'GENERAL':{'userParam01':'4.31'}})
            rnek(folder,cn,True,log_suffix="_6fR",n_procs=nps)
            ccall(('visnek fsr'+cn).split(), cwd=folder)
            
            c_pf(pf,pf,{'GENERAL':{'userParam01':'4.32'}})
            rnek(folder,cn,True,log_suffix="_6fI",n_procs=nps)
            ccall(('visnek fsi'+cn).split(), cwd=folder)

            ctoc = time.perf_counter(); cttime=ctoc-ctic
            print(f"Case finished in in {cttime:0.1f} seconds")
            print(f"                    {cttime/60:0.1f} minutes")
            print(f"                    {cttime/3600:0.2f} hours")

    toc = time.perf_counter(); ttime=toc-tic
    print(f"Script finished in in {ttime:0.2f} seconds")
    print(f"                      {ttime/60:0.1f} minutes")
    print(f"                      {ttime/3600:0.2f} hours")
