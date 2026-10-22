#!/usr/bin/env python
import sys, os, stat, re, time
from unittest import skip
from shutil import copyfile
import unittest, argparse, shutil, fileinput
from subprocess import call, PIPE, STDOUT, Popen
from subprocess import check_call as ccall 
from numpy import exp, log10

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

def Stt(Re):
    #return 0.5*(1.-exp(1.8429-0.980*log10(Re)))
    return 0.208-4.16/Re #0.198*(1.-(19.7/Re))

###############################################################
if __name__ == "__main__":
    tic = time.perf_counter()

    print('to run this script:')
    print('nohup python3 -u auton.py >>logfile 2>&1 &')
    root = os.getcwd()  # main folder loation
    print("Current working directory: {0}".format(os.getcwd()))
    base = "base"  # reference case - base case to copy
    cn = "1cyl"  # case name
    p1 = ['500'] # lopp will break when target is reached
    T = 1.0/0.125
    p2 = [T/10,T/8,T/6,T/4,T/2,T,2*T]
    nps = 6

    for i in range(len(p1)):
        for j in range(len(p2)):
            ctic = time.perf_counter()
            p2[j] = round(p2[j],3)
            path = p1[i] + "_" + str(p2[j])
            if os.path.isdir(path):
                print('Case exists:',path)
                continue

            ccopy(base,path)
            folder = root +'/' + path + '/'
            pf = folder + cn + ".par"

            c_pf(pf,pf,{'GENERAL':{'endTime':str(p2[j])}})
            c_pf(pf,pf,{'GENERAL':{'userParam07':str(int(p1[i]))}})
            rnek(folder,cn,True,log_suffix="_dir",n_procs=nps)

            ctoc = time.perf_counter(); cttime=ctoc-ctic
            print(f"Case finished in in {cttime:0.1f} seconds")
            print(f"                    {cttime/60:0.1f} minutes")
            print(f"                    {cttime/3600:0.2f} hours")

    toc = time.perf_counter(); ttime=toc-tic
    print(f"Script finished in in {ttime:0.2f} seconds")
    print(f"                      {ttime/60:0.1f} minutes")
    print(f"                      {ttime/3600:0.2f} hours")
