EXE=meshPre
DIR_BIN=bin/
DIR_TEST=test_Case/
FC=gfortran
FLAG= -O2 -g
EXOBJS=VARS_PUBLIC.o\
       READ_GRID.o\
       PERFORM_OVERLAP1.o\
       PERFORM_OVERLAP2.o\
       PERFORM_BLK_SPLIT.o\
       PERFORM_PNT_LEAP.o\
       FUNCTION_BLK_SPLIT.o\
       OUTPUT_GRID.o\
       GRIDCON_GEN.o\
       BCTYPE_GEN.o\
       RUN_GEN.o\
       INITDAMP_GEN.o\
       PSEXIT.o\
       WRITE2PLT.o\
       grid-split-glue.o

.PHONY: clean cleanall

${EXE}:${EXOBJS}
	${FC} ${FLAG} -o ${DIR_BIN}$@ ${EXOBJS}
# Comment the following 2 lines upon finishing test
${EXE}:${EXOBJS}
	${FC} ${FLAG} -o ${DIR_TEST}$@ ${EXOBJS}

${EXOBJS}: %.o: %.f90
	${FC} ${FLAG} -c $<

clean:
	@rm -f *.o  *.mod
cleanall:
	@rm -f ${DIR_BIN}${EXE} *.o *.mod
