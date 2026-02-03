#!/bin/bash

$HOME/MATLAB_R2020a/bin/matlab -nodisplay -nosplash -nodesktop -r "try, run('main.m'), catch me, fprintf('%s / %s\n',me.identifier,me.message), end, exit" > logfile
