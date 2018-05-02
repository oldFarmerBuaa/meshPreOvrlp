#!/bin/bash
EXE=`ls 2d-*`
mpirun -np  20 ./${EXE} &
