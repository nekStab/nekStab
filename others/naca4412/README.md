Fully automatic parametric study of the NACA4412 airfoil.

mshs/ contain the .re2 and .ma2 for different AoA.

To run:
python3 autobot_naca.py

or remotely:

nohup python3 -u autobot_naca.py >>logfile 2>&1 &
tail -f logfile
killall python3
