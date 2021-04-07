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

CBLACKBG  = '\33[40m'
CREDBG    = '\33[41m'
CGREENBG  = '\33[42m'
CYELLOWBG = '\33[43m'
CBLUEBG   = '\33[44m'
CVIOLETBG = '\33[45m'
CBEIGEBG  = '\33[46m'
CWHITEBG  = '\33[47m'

CGREY    = '\33[90m'
CRED2    = '\33[91m'
CGREEN2  = '\33[92m'
CYELLOW2 = '\33[93m'
CBLUE2   = '\33[94m'
CVIOLET2 = '\33[95m'
CBEIGE2  = '\33[96m'
CWHITE2  = '\33[97m'

CGREYBG    = '\33[100m'
CREDBG2    = '\33[101m'
CGREENBG2  = '\33[102m'
CYELLOWBG2 = '\33[103m'
CBLUEBG2   = '\33[104m'
CVIOLETBG2 = '\33[105m'
CBEIGEBG2  = '\33[106m'
CWHITEBG2  = '\33[107m'

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
        print("Modf.",CBOLD+CITALIC+key+CEND, "to",CBOLD+value+CEND,'in',infile)
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
            if section=='GENERAL':
                print("In",section,':',CGREEN+name+CEND,'set to',CURL+CGREEN+CBOLD+val+CEND,'in',infile)
            elif section=='VELOCITY':
                print("In",section,':',CBLUE+name+CEND,'set to',CURL+CBLUE+CBOLD+val+CEND,'in',infile)
            elif section=='PRESSURE': 
                print("In",section,':',CYELLOW+name+CEND,'set to',CURL+CYELLOW+CBOLD+val+CEND,'in',infile)
            elif section=='TEMPERATURE':
                print("In",section,':',CRED+name+CEND,'set to',CURL+CRED+CBOLD+val+CEND,'in',infile)
            else:
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

def Stt(Re):
    #return 0.5*(1.-exp(1.8429-0.980*log10(Re)))
    return 0.208-4.16/Re #0.198*(1.-(19.7/Re))

###############################################################
if __name__ == "__main__":
    tic = time.perf_counter()

    print('to run this script:')
    print('nohup python3 -u autobot.py >>logfile 2>&1 &')
    root = os.getcwd()  # main folder loation
    print("Current working directory: {0}".format(os.getcwd()))
    base = "base"  # reference case - base case to copy
    cn = "1cyl"  # case name


    aoa = ['00','10','20']
    re = ['']
    print(p1)
    nps = 6

    for i in range(len(p1)):
        for j in range(len(p2)):
            ctic = time.perf_counter()

            path = p1[i] + "_" + p2[j]
            
            ccopy(base,path)
            
            folder = root +'/' + path + '/'
            pf = folder + cn + ".par"
            SZ = folder + "SIZE"
            # {section: {optname: value, ...}, ...}
            print(pf)

            # PARAMS
            Re = float(p1[i])
            Pe = 0.71*Re
            St = abs(Stt(Re))
            Tau = round((1/St)/8,2)
            print('Re=',Re,' Pe=',Pe,' St=',St,' Tau=',Tau)

            # BASE FLOW
            tole = '1.E-7'; tolep = '1e-5'
            c_pf(pf,pf,{'GENERAL':{'endTime':300*str(Tau)}})
            c_pf(pf,pf,{'GENERAL':{'userParam01':'0'}})
            c_pf(pf,pf,{'GENERAL':{'userParam03':'0'}})
            #c_pf(pf,pf,{'GENERAL':{'userParam04':str(St)}})
            #c_pf(pf,pf,{'GENERAL':{'userParam05':'0.1'}})
            c_pf(pf,pf,{'GENERAL':{'variableDt':'yes'}})
            c_pf(pf,pf,{'GENERAL':{'targetCFL':'0.5'}})
            c_pf(pf,pf,{'GENERAL':{'extrapolation':'standard'}})
            c_pf(pf,pf,{'VELOCITY':{'viscosity':'-'+str(Re)}})
            c_pf(pf,pf,{'VELOCITY':{'residualtol':tole}})
            c_pf(pf,pf,{'VELOCITY':{'residualproj':'no'}})
            c_pf(pf,pf,{'PRESSURE':{'residualtol':tolep}})
            c_pf(pf,pf,{'PRESSURE':{'residualproj':'no'}})
            c_pf(pf,pf,{'TEMPERATURE':{'conductivity':'-'+str(Pe)}})
            c_pf(pf,pf,{'TEMPERATURE':{'residualtol':tole}})
            c_pf(pf,pf,{'TEMPERATURE':{'residualproj':'no'}})
            rnek(folder,cn,True,log_suffix="_1b",n_procs=nps)
            ccall(('visnek '+cn).split(), cwd=folder)

            # BASE FLOW
            tole = '1.E-9'
            cSZ(SZ,SZ,{"k_dim": "64"})
            c_pf(pf,pf,{'GENERAL':{'endTime':str(Tau)}})
            c_pf(pf,pf,{'GENERAL':{'userParam01':'1'}})
            c_pf(pf,pf,{'GENERAL':{'userParam03':'3'}})
            #c_pf(pf,pf,{'GENERAL':{'userParam04':str(St)}})
            #c_pf(pf,pf,{'GENERAL':{'userParam05':'0.1'}})
            c_pf(pf,pf,{'GENERAL':{'variableDt':'no'}})
            c_pf(pf,pf,{'GENERAL':{'targetCFL':'0.5'}})
            c_pf(pf,pf,{'GENERAL':{'extrapolation':'standard'}})
            c_pf(pf,pf,{'VELOCITY':{'viscosity':'-'+str(Re)}})
            c_pf(pf,pf,{'VELOCITY':{'residualtol':tole}})
            c_pf(pf,pf,{'VELOCITY':{'residualproj':'no'}})
            c_pf(pf,pf,{'PRESSURE':{'residualtol':tole}})
            c_pf(pf,pf,{'PRESSURE':{'residualproj':'no'}})
            c_pf(pf,pf,{'TEMPERATURE':{'conductivity':'-'+str(Pe)}})
            c_pf(pf,pf,{'TEMPERATURE':{'residualtol':tole}})
            c_pf(pf,pf,{'TEMPERATURE':{'residualproj':'yes'}})
            rnek(folder,cn,True,log_suffix="_1b",n_procs=nps)
            shutil.copyfile(folder+'NBF'+cn+'0.f00001',folder+'BF_'+cn+'0.f00001')
            ccall(('visnek BF_'+cn).split(), cwd=folder)

            # DIRECT 
            tole = tole ; tolep = '1.0E-7'
            cSZ(SZ,SZ,{"k_dim": "64"})
            c_pf(pf,pf,{'GENERAL':{'endTime':str(Tau)}})
            c_pf(pf,pf,{'GENERAL':{'userParam01':'3'}})
            c_pf(pf,pf,{'GENERAL':{'userParam03':'0'}})
            c_pf(pf,pf,{'GENERAL':{'userParam08':'5'}})
            c_pf(pf,pf,{'GENERAL':{'userParam09':'5'}})
            c_pf(pf,pf,{'GENERAL':{'userParam10':'1'}})
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
            # O is changed to v in code -- usrdat2
            c_pf(pf,pf,{'GENERAL':{'userParam08':'5'}})
            c_pf(pf,pf,{'GENERAL':{'userParam09':'5'}})
            c_pf(pf,pf,{'GENERAL':{'userParam10':'1.7'}})
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

            ctoc = time.perf_counter(); cttime=ctoc-ctic
            print(f"Case finished in in {cttime:0.1f} seconds")
            print(f"                    {cttime/60:0.1f} minutes")
            print(f"                    {cttime/3600:0.2f} hours")

    toc = time.perf_counter(); ttime=toc-tic
    print(CBLINK + f"Script finished in in {ttime:0.2f} seconds" + CEND)
    print(CBLINK + f"                      {ttime/60:0.1f} minutes" + CEND)
    print(CBLINK + f"                      {ttime/3600:0.2f} hours" + CEND)
