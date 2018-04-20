#!/bin/bash
INITIAL_SLEEP=10
EXE=`ls 2d-*`
mpirun -np  35 -x ${INITIAL_SLEEP} ./${EXE} &
