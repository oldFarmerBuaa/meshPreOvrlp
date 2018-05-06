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
INTEGER i_BlkGlb,i_Point,i_Blk,max_nPntFnl,min_LvOffUnfrmGlb,max_nPntEqvFnl
INTEGER i_Proc,n_PntByNow,flag_SortByBubble,ratio_Balance
INTEGER n_Loop,is_StepUp,idx_GroupBase
REAL(8) width_GroupInit,width_StepInit,width_Group,width_Step,diff_Ratio,dt_Base
INTEGER,DIMENSION(:),ALLOCATABLE :: min_LvOffUnfrmFnl,idx_Bubble
INTEGER,DIMENSION(:),ALLOCATABLE :: min_LvOffUnfrmFnlLow2High
INTEGER,DIMENSION(:),ALLOCATABLE :: var_In,n_PntFinalLow2High,Lv_Recorded,n_PntEqvFnlLow2High
REAL(8),DIMENSION(:),ALLOCATABLE :: var_InReal,dt_FinalLow2High,dt_Low2HighNormalized
ALLOCATE(II_Final(n_BlkFinal),JJ_Final(n_BlkFinal),KK_Final(n_BlkFinal))
ALLOCATE(n_PntFinal(n_BlkFinal),idx_Proc(n_BlkFinal),idx_Mltstp(n_BlkFinal))
ALLOCATE(n_PntEqvFnl(n_BlkFinal),dt_Final(n_BlkFinal))
ALLOCATE(Lv_OffUnfrmFnl(3,n_BlkFinal))
ALLOCATE(min_LvOffUnfrmFnl(n_BlkFinal))
ALLOCATE(min_LvOffUnfrmFnlLow2High(n_BlkFinal))
ALLOCATE(var_In(n_BlkFinal),idx_Bubble(n_BlkFinal))
ALLOCATE(var_InReal(n_BlkFinal))
ALLOCATE(dt_FinalLow2High(n_BlkFinal))
ALLOCATE(dt_Low2HighNormalized(n_BlkFinal))
ALLOCATE(n_PntFinalLow2High(n_BlkFinal))
ALLOCATE(n_PntEqvFnlLow2High(n_BlkFinal))
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
  READ(30,*) dt_Final(J),n_PntEqvFnl(J)
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
max_nPntFnl       = MAXVAL(n_PntFinal)
max_nPntEqvFnl    = MAXVAL(n_PntEqvFnl)
n_PntFnlTot       = SUM(n_PntFinal)
n_PntEqvFnlTot    = SUM(n_PntEqvFnl)
min_LvOffUnfrmGlb = MINVAL(min_LvOffUnfrmFnl)

!---------------------------------------------------------------------------!
! Calculate No. time steps
!var_In=min_LvOffUnfrmFnl
!CALL bubble_Int(var_In,min_LvOffUnfrmFnlLow2High,idx_Bubble,n_BlkFinal)
!Lv_Recorded=100
!Lv_Mltstp=0
!DO J=1,n_BlkFinal
!  tmp_Int=0
!  DO J2=1,n_BlkFinal
!    IF(min_LvOffUnfrmFnlLow2High(J) .NE. Lv_Recorded(J2)) THEN
!      tmp_Int=tmp_Int+1
!    ENDIF
!  ENDDO
!  ! Not in record
!  IF(tmp_Int==n_BlkFinal) THEN
!    Lv_Mltstp=Lv_Mltstp+1
!    ! Put in record
!    Lv_Recorded(J)=min_LvOffUnfrmFnlLow2High(J)
!  ENDIF
!ENDDO

!---------------------------------------------------------------------------!
! Calculate idx_Mltstp
!I=0
!tmp_Int=100
!DO J=1,n_BlkFinal
!   IF(min_LvOffUnfrmFnlLow2High(J) .NE. tmp_Int) THEN
!     I=I+1
!     tmp_Int=min_LvOffUnfrmFnlLow2High(J)
!   ENDIF   
!   idx_Mltstp(idx_Bubble(J))=I
!ENDDO

!---------------------------------------------------------------------------!
! Calculate Lv_OffUnfrmFnl for all blocks after split
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

!---------------------------------------------------------------------------!
! Group time steps: dt_Final and Assign idx_Mltstp
! The time step index assignment is linear by now
var_InReal=dt_Final
CALL bubble_Real(var_InReal,dt_FinalLow2High,idx_Bubble,n_BlkFinal)
WRITE(*,*) "dt_Final:"
WRITE(*,200) dt_Final
WRITE(*,*) "dt_FinalLow2High:"
WRITE(*,200) dt_FinalLow2High
200 FORMAT(4F20.15)
WRITE(*,*) "Lv_Mltstp:",Lv_Mltstp
idx_Mltstp=0
width_GroupInit=0.05
width_StepInit=0.025
width_Group=width_GroupInit
width_Step=width_StepInit
idx_Mltstp(idx_Bubble(n_BlkFinal))=0
n_Loop=0
is_StepUp=1
DO WHILE(idx_Mltstp(idx_Bubble(n_BlkFinal)) .NE. Lv_Mltstp)
  ! Upon finish one cycle check if max_idx_Mltstp = Lv_Mltstp
  ! If max_idx_Mltstp > Lv_Mltstp, increase width_Group, otherwise reduce it and re-run
  n_Loop=n_Loop+1
  IF(n_Loop >1 .AND. idx_Mltstp(idx_Bubble(n_BlkFinal)) > Lv_Mltstp) THEN
    IF(is_StepUp==0) width_Step=width_Step/2
    is_StepUp=1
    width_Group=width_Group+width_Step
  ELSEIF(n_Loop > 1 .AND. idx_Mltstp(idx_Bubble(n_BlkFinal)) < Lv_Mltstp) THEN
    is_StepUp=0
    DO WHILE(width_Group < width_Step)
      width_Step=width_Step/2
    ENDDO
    width_Group = width_Group-width_Step
  ENDIF
  DO J=1,n_BlkFinal
    dt_Low2HighNormalized(J)=dt_FinalLow2High(J)/dt_FinalLow2High(n_BlkFinal)
    IF(J==1) THEN 
      idx_Mltstp(idx_Bubble(J))=1
      dt_Base       = dt_Low2HighNormalized(J)
      idx_GroupBase = idx_Mltstp(idx_Bubble(J))
    ELSE
      diff_Ratio=(dt_Low2HighNormalized(J)-dt_Base)/dt_Base
      IF(diff_Ratio > width_Group) THEN
        ! If diff_Ratio > width_Group, set group index to be idx_Group_Base+1
        ! And set this block as base block
        idx_Mltstp(idx_Bubble(J))=idx_GroupBase+1
        idx_GroupBase = idx_Mltstp(idx_Bubble(J))
        dt_Base       = dt_Low2HighNormalized(J)
      ELSE    
        ! If diff_Ratio <= width_Group, set the same group index with base
        idx_Mltstp(idx_Bubble(J))=idx_GroupBase
      ENDIF
    ENDIF
  ENDDO
ENDDO
WRITE(*,*) "dt_Low2HighNormalized:"
WRITE(*,200) dt_Low2HighNormalized
WRITE(*,*) "idx_Mltstp::"
WRITE(*,300) idx_Mltstp
300 FORMAT(4I2)
!---------------------------------------------------------------------------!
! Assign processors with ratio_Blance give as below
! ratio_Balance=1.5 means the max block can be 1.5 times larger than the min.
!!! Memory balance is not garanteed if bubble up equivalent points
!!! It's better to apply Jiangming's strategy to achieve CPU-time and Memory balance simultaneously
ratio_Balance=1.5
IF(flag_Mltstp==1)  THEN
  ! Assign based on equivalent points
  n_PntByNow=0
  i_Proc=1
  IF(flag_SortByBubble==1) THEN
    ! ASSIGN IDX_PROC SORTING BY POINTS
    var_In=n_PntEqvFnl
    CALL bubble_Int(var_In,n_PntEqvFnlLow2High,idx_Bubble,n_BlkFinal)
    DO J=1,n_BlkFinal
      n_PntByNow=n_PntByNow+n_PntEqvFnlLow2High(J)
      IF(n_PntByNow<ratio_Balance*max_nPntEqvFnl) THEN
        idx_Proc(idx_Bubble(J))=i_Proc
      ELSE
        i_Proc=i_Proc+1
        idx_Proc(idx_Bubble(J))=i_Proc
        n_PntByNow=n_PntEqvFnlLow2High(J)
      ENDIF
    ENDDO
  ELSE
    ! ASSIGN IDX_PROC WITHOUT SORTING BY POINTS
    DO J=1,n_BlkFinal
    !!! Sort Block index by n_PntFinal will be better
      n_PntByNow=n_PntByNow+n_PntEqvFnl(J)
      IF(n_PntByNow<ratio_Balance*max_nPntEqvFnl) THEN
        idx_Proc(J)=i_Proc
      ELSE
        i_Proc=i_Proc+1
        idx_Proc(J)=i_Proc
        n_PntByNow=n_PntEqvFnl(J)
      ENDIF
    ENDDO
  ENDIF
ELSE
  ! Assign based on physicial/real points
  n_PntByNow=0
  i_Proc=1
  IF(flag_SortByBubble==1) THEN
    ! ASSIGN IDX_PROC SORTING BY POINTS
    var_In=n_PntFinal
    CALL bubble_Int(var_In,n_PntFinalLow2High,idx_Bubble,n_BlkFinal)
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
                  n_PntFinal(J),idx_Mltstp(J),n_PntEqvFnl(J)
  ELSE
    WRITE(30,100) J,II_Final(J),JJ_Final(J),KK_Final(J),&
    idx_Proc(J),n_PntFinal(J),idx_Mltstp(J),n_PntEqvFnl(J)
  ENDIF
ENDDO
WRITE(30,*) "    IMB	I	J	idx_Proc n_PntFinal idx_Mltstp n_PntEqvFnl"
WRITE(30,*) ""
WRITE(30,*) "             Total points: ",n_PntFnlTot
WRITE(30,*) "  Total equivalent points: ",n_PntEqvFnlTot
WRITE(30,*) "Total processors assigned: ",max_Proc
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
                  n_PntFinal(J),J,n_PntEqvFnl(J)
  ELSE
    WRITE(30,100) J,II_Final(J),JJ_Final(J),KK_Final(J),&
    idx_Proc(J),n_PntFinal(J),J,n_PntEqvFnl(J)
  ENDIF
ENDDO
WRITE(30,*) "    IMB	I	J	idx_Proc n_PntFinal idx_Mltstp n_PntEqvFnl"
WRITE(30,*) ""
WRITE(30,*) "             Total points: ",n_PntFnlTot
WRITE(30,*) "  Total equivalent points: ",n_PntEqvFnlTot
WRITE(30,*) "Total processors assigned:",max_Proc
CLOSE(30)
100 FORMAT(3X,I4,3X,I6,3X,I6,3X,I3,3X,I8,3X,I3,3X,I6)
END


!---------------------------------------------------------------------------!
SUBROUTINE bubble_Int(var_In,var_Out,idx_Bubble,n_Dim)
! Bubble up sorting
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


!---------------------------------------------------------------------------!
SUBROUTINE bubble_Real(var_In,var_Out,idx_Bubble,n_Dim)
! Bubble up sorting
IMPLICIT NONE

INTEGER n_Dim,I,J
INTEGER idx_Bubble(n_Dim),idx_BubbleTmp(n_Dim)
REAL(8) tmp1,tmp2,var_In(n_Dim),var_Out(n_Dim)
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
