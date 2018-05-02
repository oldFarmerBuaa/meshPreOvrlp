#!/bin/bash
export INITIAL_SLEEP_TIME=10
EXE=`ls 2d-*`
mpirun -np  20 -x INITIAL_SLEEP_TIME ./${EXE} &
