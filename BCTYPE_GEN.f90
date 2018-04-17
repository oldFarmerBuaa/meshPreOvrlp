subroutine bctype_Gen
!
USE read_gridcon_pre
USE read_mesh_original
USE system_var
USE point_leaping
USE block_overlaping
USE block_splitting
!
IMPLICIT NONE
INTEGER I,J
REAL(8) xc,yc
ALLOCATE(isidtp_Final(6,n_BlkFinal))

xc=0.0
yc=0.0
101 FORMAT(7X,F3.1,7X,F3.1,7X,I4,7X,I4,7X,I4,7X,I4,7X,I3)
102 FORMAT(7X,F3.1,7X,F3.1,7X,I4,7X,I4,7X,I4,7X,I4,7X,I4,7X,I4,7X,I3)
!-------------------------------------------------------------------------------!
DO J=1,n_BlkFinal
  WRITE(FILENAME,'(A10,I4.4)') "GRIDBLO_BC",J
  OPEN(30,FILE=FILENAME,STATUS="OLD")
  READ(30,*)
  READ(30,*)
  READ(30,*)
  READ(30,*) isidtp_Final(1,J),isidtp_Final(2,J),isidtp_Final(3,J),&
             isidtp_Final(4,J),isidtp_Final(5,J),isidtp_Final(6,J)
  CLOSE(30)
ENDDO
!-------------------------------------------------------------------------------!
WRITE(FILENAME,'(A18)') "output_GRID/BCTYPE"
OPEN(30,FILE=FILENAME,STATUS="UNKNOWN")
IF(cas_Dim=="2D") THEN
 ! WRITE(30,'(3X,A6,3X,A3,3X,A3,3X,A3,3X,A3,3X,A8)') "ISIDTP",&
 !          "(1)","(2)","(3)","(4)","idx_Proc"
  WRITE(30,*) n_BlkFinal
  DO J=1,n_BlkFinal
     WRITE(30,*) J,(isidtp_Final(I,J),I=1,4),idx_Proc(J)
  ENDDO
ELSE
 ! WRITE(30,'(3X,A6,3X,A3,3X,A3,3X,A3,3X,A3,3X,A3,3X,A3,3X,A8)') "ISIDTP",&
 !          "(1)","(2)","(3)","(4)","(5)","(6)","idx_Proc"
  WRITE(30,*) n_BlkFinal
  DO J=1,n_BlkFinal
     WRITE(30,*) J,(isidtp_Final(I,J),I=1,6),idx_Proc(J)
  ENDDO
ENDIF  
CLOSE(30)
!-------------------------------------------------------------------------------!
WRITE(FILENAME,'(A22)') "output_GRID/BC4CASECON"
OPEN(30,FILE=FILENAME,STATUS="UNKNOWN")
IF(cas_Dim=="2D") THEN
  DO J=1,n_BlkFinal
     WRITE(30,101) xc,yc,(isidtp_Final(I,J),I=1,4),idx_Proc(J)
  ENDDO
ELSE
  DO J=1,n_BlkFinal
     WRITE(30,102) xc,yc,(isidtp_Final(I,J),I=1,6),idx_Proc(J)
  ENDDO
ENDIF  
CLOSE(30)
END
