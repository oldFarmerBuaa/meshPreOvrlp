#!/bin/bash
DIR_BIN=bin/
gfortran -g -o ${DIR_BIN}wall_bcGen wall_bcGen.f90
gfortran -g -o ${DIR_BIN}write2plt_input write2plt_input.f90
gfortran -g -o ${DIR_BIN}write2plt_ovrlp write2plt_ovrlp.f90
gfortran -g -o ${DIR_BIN}write2plt_pntlp write2plt_pntlp.f90
gfortran -g -o ${DIR_BIN}write2plt_final write2plt_final.f90
