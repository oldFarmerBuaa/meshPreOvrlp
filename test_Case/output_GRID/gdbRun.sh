#!/bin/bash
export INITIAL_SLEEP_TIME=10
EXE=`ls 2d-*`
mpirun -np  35 -x INITIAL_SLEEP_TIME ./${EXE} &
