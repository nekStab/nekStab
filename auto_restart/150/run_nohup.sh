#!/bin/bash
# Run run.sh in the background with nohup. The output of run.sh goes to
# logfile, logerror and logrestart, so nohup.out would stay empty: discard it.

nohup ./run.sh > /dev/null 2>&1 &
