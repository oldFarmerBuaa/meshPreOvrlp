#!/bin/bash
EXE=`ls 2d-*`
mpirun -np  21 ./${EXE} &
