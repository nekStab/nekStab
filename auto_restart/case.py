import sys
sys.path.append("/gpfswork/rech/vpo/rvpo014")
from autotools import *

cn = "1cyl"

pbs_file = 'jz.pbs'
next_script = 'check_next'
file_ct = next_script + '.cc'
file_ft = next_script + '.ft'
cns = cn[0]
bf = "BF_"+cn+"0.f00001"
liftdrag = "lift_drag.dat"
hisfile = cn+".his"
globenergy = "total_energy.dat"
globenstro = "total_enstrophy.dat"
final_dns_file = cn+"0.f00002"

tol = 1e-10
def adjust_par_file(pf, Re):
    """Adjust parameters in the .par file for the simulation."""
    params = [
        ('GENERAL', 'startFrom', bf),
        ('GENERAL', 'endTime', '1.25'),
        ('GENERAL', 'userParam01', '2'),
        ('GENERAL', 'userParam07', '200'),
        ('GENERAL', 'userParam10', '1.7'),
        ('PRESSURE', 'residualtol', str(tol)),
        ('VELOCITY', 'residualtol', str(tol)),
        ('VELOCITY', 'viscosity', f"-{float(Re)}")
    ]
    for section, key, value in params:
        c_pf(pf, pf, {section: {key: value}})
