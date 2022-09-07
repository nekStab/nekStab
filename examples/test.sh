#/bin/bash
nohup python3 -u test.py > logfile 2>&1 & # overwrite
#nohup python3 -u test.py >> logfile 2>&1 & # append
