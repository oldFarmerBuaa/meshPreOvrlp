subroutine run_Gen
!
USE system_var,ONLY: max_Proc,FILENAME
!
IMPLICIT NONE
! Case run file
WRITE(FILENAME,'(A21)') "output_GRID/casRun.sh"
OPEN(30,FILE=FILENAME,STATUS="UNKNOWN")
WRITE(30,'(A11)') "#!/bin/bash"
WRITE(30,'(A13)') "EXE=`ls 2d-*`"
WRITE(30,'(A11,I3,A11)') "mpirun -np ",max_Proc," ./${EXE} &"
CLOSE(30)
CALL system('chmod 755 output_GRID/casRun.sh')
! Debug run file
WRITE(FILENAME,'(A21)') "output_GRID/dbgRun.sh"
OPEN(30,FILE=FILENAME,STATUS="UNKNOWN")
WRITE(30,'(A11)') "#!/bin/bash"
WRITE(30,'(A16)') "INITIAL_SLEEP=10"
WRITE(30,'(A13)') "EXE=`ls 2d-*`"
WRITE(30,'(A11,I3,A31)') "mpirun -np ",max_Proc," -x ${INITIAL_SLEEP} ./${EXE} &"
CLOSE(30)
CALL system('chmod 755 output_GRID/dbgRun.sh')
RETURN
END SUBROUTINE
