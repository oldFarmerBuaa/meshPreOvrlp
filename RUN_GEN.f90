subroutine run_Gen
!
USE system_var,ONLY: max_Proc,FILENAME
!
IMPLICIT NONE
WRITE(FILENAME,'(A23)') "output_GRID/run_Case.sh"
OPEN(30,FILE=FILENAME,STATUS="UNKNOWN")
WRITE(30,'(A11)') "#!/bin/bash"
WRITE(30,'(A13)') "EXE=`ls 2d-*`"
WRITE(30,'(A10,I3,A11)') "mpirun -n ",max_Proc," ./${EXE} &"
CALL system('chmod 755 output_GRID/run_Case.sh')
RETURN
END SUBROUTINE
