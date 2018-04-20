#!/bin/bash
EXE=`ls 2d-*`
mpirun -np  35 ./${EXE} &
