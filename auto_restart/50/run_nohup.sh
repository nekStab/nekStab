#!/bin/bash
# Script to run run.pbs in background with nohup, sending all output to logfile

nohup ./run.sh & # > logfile 2> logerror &