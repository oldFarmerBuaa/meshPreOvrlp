subroutine initdamp_Gen
!---------------------------------------------------------------!
! Generate initial background damping coef., b.c. damping coef.	!
! and newtonian damping						!
!								!
!---------------------------------------------------------------!
USE read_gridcon_pre
USE read_mesh_original
USE system_var
USE point_leaping
USE block_overlaping
USE block_splitting
!
IMPLICIT NONE
INTEGER,PARAMETER :: NSPPT=2
INTEGER I,J,IMB
INTEGER NOUTBC
REAL(8) ADAMP_RE,ADAMPSB,COENT,ADAMPPML,BETA,AX,AY
REAL(8) ratio_EnhancePnt
REAL(8) XSPPT(NSPPT),YSPPT(NSPPT),ZSPPT(NSPPT)
INTEGER,DIMENSION(:),ALLOCATABLE   :: IOUTBCCON
REAL(8),DIMENSION(:),ALLOCATABLE   :: WTHSPPT,ADAMPSPPT,COENTSPPT
REAL(8),DIMENSION(:),ALLOCATABLE   :: WTHOUTBC,ADAMPOUTBC,COENTOUTBC
REAL(8),DIMENSION(:),ALLOCATABLE   :: XOUTBC,YOUTBC
REAL(8),DIMENSION(:,:),ALLOCATABLE :: ADAMPBC_RE,COENTBC,ratio_Enhance
! 4 FOR 2D CASE, 6 FOR 3D CASE
ALLOCATE(ADAMPBC_RE(n_BlkFinal,6),COENTBC(n_BlkFinal,6))
ALLOCATE(ratio_Enhance(n_BlkFinal,6))
!---------------------------------------------------------------!
! GENERATE 2D INITIAL DAMPING COEFFICIENTS
ADAMP_RE = 4.D-2
ADAMPSB  = 0.D0
COENT    = 0.D0
ADAMPPML = 2.D0
BETA     = .22331554941283124871
AX       = 0.D0
AY       = 0.D0
!---------------------------------------------------------------!
! CALCULATE DAMP ALONG BOUNDARIES
DO J=1,n_BlkFinal
  DO I=1,6
    ! MORE DAMP ON WALL, LESS DAMP ON FLUID BOUNDARY
    ! WALL FAR FROM CORE-WALLS SHOULD HAVE LESS DAMP
    IF(isidtp_Final(I,J)==7000) THEN
      ! FROM CORE-WALL TO FAR-WALL, THIS RATIO RANGE FROM 0.5D0 TO 2.D0
      ratio_Enhance(J,I)=2.D0  
    ELSE
      ratio_Enhance(J,I)=5.D-1
    ENDIF
    ADAMPBC_RE(J,I)=ratio_Enhance(J,I)*ADAMP_RE
    COENTBC(J,I)=ratio_Enhance(J,I)*COENT
  ENDDO
ENDDO
!---------------------------------------------------------------!
! ALLOCATE ARRAY FOR EXTRA DAMP ALONG GIVING LINE
NOUTBC = 0
ALLOCATE(IOUTBCCON(NOUTBC))
ALLOCATE(XOUTBC(NOUTBC), YOUTBC(NOUTBC),&
         WTHOUTBC(NOUTBC), ADAMPOUTBC(NOUTBC), COENTOUTBC(NOUTBC))
DO I=1,NOUTBC
   ! SET VALUE FOR EXTRA LINES, WHAT'S THE RULE?
ENDDO
!---------------------------------------------------------------!
! CALCULATE DAMP AT EXTRA POINTS
ALLOCATE(WTHSPPT(NSPPT),ADAMPSPPT(NSPPT),COENTSPPT(NSPPT))
DATA XSPPT /-0.5, 0.5/
DATA YSPPT / 0.5, 0.5/
DATA ZSPPT / 0.0, 0.0/
WTHSPPT=5.D-2
ratio_EnhancePnt=2.D0
ADAMPSPPT=ratio_EnhancePnt*ADAMP_RE
COENTSPPT=ratio_EnhancePnt*COENT
!---------------------------------------------------------------!
! WRITE TO FILE INIDAMP
OPEN(10,FILE='output_GRID/INIDAMP',STATUS='UNKNOWN')
WRITE(10,*) "ADAMP_RE   ADAMPSB   COENT   ADAMPPML   BETA   AX   AY"
! Write background damping coef.
WRITE(10,50) ADAMP_RE, ADAMPSB, COENT, ADAMPPML, BETA, AX, AY
! Write damping coef. along every b.c.
WRITE(10,200) "ADAMPBC_RE","(1)","(2)","(3)","(4)","IMB"
DO IMB = 1,n_BlkFinal
  WRITE(10,100) ADAMPBC_RE(IMB,1),ADAMPBC_RE(IMB,2),&
             ADAMPBC_RE(IMB,3),ADAMPBC_RE(IMB,4),IMB
ENDDO
! Write  newtonian cooling type damping coef. along every b.c.
WRITE(10,200) "COENTBC","(1)","(2)","(3)","(4)","IMB"
DO IMB = 1,n_BlkFinal
  WRITE(10,100) COENTBC(IMB,1),COENTBC(IMB,2),&
                COENTBC(IMB,3),COENTBC(IMB,4),IMB
ENDDO
! Write extra added damping along those lines.
WRITE(10,*) "NUMBER OF EXTRA DAMPING ALONG LINES: NOUTBC"
WRITE(10,*) NOUTBC
WRITE(10,*) "IOUTBCCON   XOUTBC    YOUTBC   WTHOUTBC    ADAMPOUTBC   COENTOUTBC" 
DO I= 1,NOUTBC
  WRITE(10,*)IOUTBCCON(I),XOUTBC(I),YOUTBC(I),&
             WTHOUTBC(I),ADAMPOUTBC(I),COENTOUTBC(I)
ENDDO
! Write extra added damping in points.
WRITE(10,*) "NSPPT"
WRITE(10,*) NSPPT 
WRITE(10,*) "XSPPT    YSPPT    WITHSPPT    ADAMPSPPT    COENTSPPT"
DO I = 1,NSPPT
  WRITE(10,150) XSPPT(I),YSPPT(I),WTHSPPT(I),ADAMPSPPT(I),COENTSPPT(I)
ENDDO
CLOSE(10)
50  FORMAT(D8.1,2X,D8.1,2X,D8.1,2X,D8.1,2X,F22.20,2X,D8.1,2X,D8.1)
100 FORMAT(12X,D8.1,2X,D8.1,2X,D8.1,2X,D8.1,2X,I3)
150 FORMAT(D8.1,2X,D8.1,2X,D8.1,2X,D8.1,2X,2X,D8.1)
200 FORMAT(A10,4X,A3,7X,A3,7X,A3,7X,A3,6X,A3)
END
