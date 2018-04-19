#!/bin/bash
EXE=`ls 2d-*`
mpirun -n  35 ./${EXE} &
