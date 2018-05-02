subroutine gridcon_Gen
!
USE read_gridcon_pre
USE read_mesh_original
USE system_var
USE point_leaping
USE block_overlaping
USE block_splitting
!
IMPLICIT NONE
INTEGER I,J,K,J2,tmp
INTEGER i_BlkGlb,i_Point,i_Blk,max_nPntFnl,min_LvOffUnfrmGlb
INTEGER i_Proc,n_PntByNow,flag_SortByBubble,ratio_Balance
INTEGER,DIMENSION(:),ALLOCATABLE :: min_LvOffUnfrmFnl,idx_Bubble
INTEGER,DIMENSION(:),ALLOCATABLE :: min_LvOffUnfrmFnlLow2High
INTEGER,DIMENSION(:),ALLOCATABLE :: var_In,n_PntFinalLow2High,Lv_Recorded
ALLOCATE(II_Final(n_BlkFinal),JJ_Final(n_BlkFinal),KK_Final(n_BlkFinal))
ALLOCATE(n_PntFinal(n_BlkFinal),idx_Proc(n_BlkFinal),idx_Mltstp(n_BlkFinal))
ALLOCATE(Lv_OffUnfrmFnl(3,n_BlkFinal))
ALLOCATE(min_LvOffUnfrmFnl(n_BlkFinal))
ALLOCATE(min_LvOffUnfrmFnlLow2High(n_BlkFinal))
ALLOCATE(var_In(n_BlkFinal),idx_Bubble(n_BlkFinal))
ALLOCATE(n_PntFinalLow2High(n_BlkFinal))
ALLOCATE(Lv_Recorded(n_BlkFinal))

! Sort points on each block by bubble up method(0/1=no/yes)
flag_SortByBubble=1

DO J=1,n_BlkFinal
  WRITE(FILENAME,'(A10,I4.4)') "GRIDBLO_BC",J
  OPEN(30,FILE=FILENAME,STATUS="OLD")
  READ(30,*)
  READ(30,*) (Lv_OffUnfrmFnl(I,J),I=1,3)
  READ(30,*)
  READ(30,*)
  READ(30,*)
  IF(cas_Dim=="2D") THEN
    READ(30,*) JJ_Final(J),II_Final(J)
    n_PntFinal(J)=II_Final(J)*JJ_Final(J)
  ELSE
    READ(30,*) JJ_Final(J),JJ_Final(J),II_Final(J)
    n_PntFinal(J)=II_Final(J)*JJ_Final(J)*KK_Final(J)
  ENDIF
  CLOSE(30)
ENDDO
DO J=1,n_BlkFinal
  tmp=100
  IF(cas_Dim=="2D") THEN
    DO I=1,2
      IF(Lv_OffUnfrmFnl(I,J)<tmp) THEN
        tmp=Lv_OffUnfrmFnl(I,J)
      ENDIF
    ENDDO
  ELSE
    DO I=1,3
      IF(Lv_OffUnfrmFnl(I,J)<tmp) THEN
        tmp=Lv_OffUnfrmFnl(I,J)
      ENDIF
    ENDDO
  ENDIF
  min_LvOffUnfrmFnl(J)=tmp
ENDDO
max_nPntFnl=MAXVAL(n_PntFinal)
min_LvOffUnfrmGlb=MINVAL(min_LvOffUnfrmFnl)
!---------------------------------------------------------------------------!
var_In=min_LvOffUnfrmFnl
CALL bubble_Int(var_In,min_LvOffUnfrmFnlLow2High,idx_Bubble,n_BlkFinal)
! Calculate No. time steps
Lv_Recorded=100
Lv_Mltstp=0
DO J=1,n_BlkFinal
  tmp_Int=0
  DO J2=1,n_BlkFinal
    IF(min_LvOffUnfrmFnlLow2High(J) .NE. Lv_Recorded(J2)) THEN
      tmp_Int=tmp_Int+1
    ENDIF
  ENDDO
  ! Not in record
  IF(tmp_Int==n_BlkFinal) THEN
    Lv_Mltstp=Lv_Mltstp+1
    ! Put in record
    Lv_Recorded(J)=min_LvOffUnfrmFnlLow2High(J)
  ENDIF
ENDDO
!---------------------------------------------------------------------------!
! Calculate idx_Mltstp
I=0
tmp_Int=100
DO J=1,n_BlkFinal
   IF(min_LvOffUnfrmFnlLow2High(J) .NE. tmp_Int) THEN
     I=I+1
     tmp_Int=min_LvOffUnfrmFnlLow2High(J)
   ENDIF   
   idx_Mltstp(idx_Bubble(J))=I
ENDDO
idx_Proc=0
i_Point=0
DO K=1,n_BlkSupPre
  DO J=1,n_BlkPre(K)
    DO I=1,2**n_Split(J,K)
      i_Blk=i_Point+I
      Lv_OffUnfrmFnl(1:3,i_Blk)=Lv_OffUniform(1:3,J,K) 
    ENDDO
    i_Point=i_Blk
  ENDDO
ENDDO
!WRITE(*,'(3X,I2,3X,I2,3X,I2)') ((Lv_OffUnfrmFnl(I,J),I=1,3),J=1,n_BlkFinal)
!---------------------------------------------------------------------------!
n_PntFnlTot=0
DO J=1,n_BlkFinal
  n_PntFnlTot=n_PntFnlTot+n_PntFinal(J)
ENDDO
n_PntByNow=0
i_Proc=1
IF(flag_SortByBubble==1) THEN
  ! ASSIGN IDX_PROC SORTING BY POINTS
  var_In=n_PntFinal
  CALL bubble_Int(var_In,n_PntFinalLow2High,idx_Bubble,n_BlkFinal)
  ratio_Balance=1.5
  DO J=1,n_BlkFinal
    n_PntByNow=n_PntByNow+n_PntFinalLow2High(J)
    IF(n_PntByNow<ratio_Balance*max_nPntFnl) THEN
      idx_Proc(idx_Bubble(J))=i_Proc
    ELSE
      i_Proc=i_Proc+1
      idx_Proc(idx_Bubble(J))=i_Proc
      n_PntByNow=n_PntFinalLow2High(J)
    ENDIF
  ENDDO
ELSE
  ! ASSIGN IDX_PROC WITHOUT SORTING BY POINTS
  DO J=1,n_BlkFinal
  !!! Sort Block index by n_PntFinal will be better
    n_PntByNow=n_PntByNow+n_PntFinal(J)
    IF(n_PntByNow<ratio_Balance*max_nPntFnl) THEN
      idx_Proc(J)=i_Proc
    ELSE
      i_Proc=i_Proc+1
      idx_Proc(J)=i_Proc
      n_PntByNow=n_PntFinal(J)
    ENDIF
  ENDDO
ENDIF
max_Proc=i_Proc
!---------------------------------------------------------------------------!
! Write GRIDCON file for calculation
WRITE(FILENAME,'(A19)') "output_GRID/GRIDCON"
OPEN(30,FILE=FILENAME,STATUS="UNKNOWN")
WRITE(30,*) n_BlkFinal,Lv_Mltstp
DO J=1,n_BlkFinal
  IF(cas_Dim=="2D") THEN
    WRITE(30,100) J,II_Final(J),JJ_Final(J),idx_Proc(J),&
                  n_PntFinal(J),idx_Mltstp(J)
  ELSE
    WRITE(30,100) J,II_Final(J),JJ_Final(J),KK_Final(J),&
    idx_Proc(J),n_PntFinal(J),idx_Mltstp(J)
  ENDIF
ENDDO
WRITE(30,*) "Total points:",n_PntFnlTot
WRITE(30,*) "Total processors assigned:",max_Proc
CLOSE(30)
!---------------------------------------------------------------------------!
! Write GRIDCON file for pre-processing to get exchange table
! Set all blocks have different time step scale
WRITE(FILENAME,'(A27)') "output_GRID/GRIDCON_EXCCTAB"
OPEN(30,FILE=FILENAME,STATUS="UNKNOWN")
WRITE(30,*) n_BlkFinal,n_BlkFinal
DO J=1,n_BlkFinal
  IF(cas_Dim=="2D") THEN
    WRITE(30,100) J,II_Final(J),JJ_Final(J),idx_Proc(J),&
                  n_PntFinal(J),J
  ELSE
    WRITE(30,100) J,II_Final(J),JJ_Final(J),KK_Final(J),&
    idx_Proc(J),n_PntFinal(J),J
  ENDIF
ENDDO
WRITE(30,*) "Total points:",n_PntFnlTot
WRITE(30,*) "Total processors assigned:",max_Proc
CLOSE(30)
100 FORMAT(3X,I4,3X,I6,3X,I6,3X,I3,3X,I8,3X,I3)
END
!---------------------------------------------------------------------------!
SUBROUTINE bubble_Int(var_In,var_Out,idx_Bubble,n_Dim)
IMPLICIT NONE

INTEGER n_Dim,I,J,tmp1,tmp2
INTEGER var_In(n_Dim),var_Out(n_Dim),idx_Bubble(n_Dim),idx_BubbleTmp(n_Dim)
! var_In = High to Low, var_Out = Low to High
DO I=1,n_Dim
  idx_Bubble(I)=I
ENDDO
DO I=1,n_Dim-1
  DO J=I+1,n_Dim
    IF(var_In(I)<var_In(J)) THEN
      tmp1=var_In(I)
      var_In(I)=var_In(J)
      var_In(J)=tmp1
      tmp2=idx_Bubble(I)
      idx_Bubble(I)=idx_Bubble(J)
      idx_Bubble(J)=tmp2
    ENDIF
  ENDDO
ENDDO
idx_BubbleTmp=idx_Bubble
DO I=1,n_Dim
  var_Out(I)=var_In(n_Dim-I+1)
  idx_Bubble(I)=idx_BubbleTmp(n_Dim-I+1)
ENDDO
RETURN
END SUBROUTINE
