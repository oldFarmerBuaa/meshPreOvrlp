program wall_bcGen
!-------------------------------------------------------------------------------!
!	Generate wall boundary parameters form original mesh generation		!
!										!
!				Xiaolong Tang 2018/03/16			!
!-------------------------------------------------------------------------------!

IMPLICIT NONE

INTEGER I,J,K
INTEGER tmp_Int
INTEGER n_LvDown,n_BlkLyr,n_LvTot,ratio_ynVSy1
INTEGER n_PntOvrlp,n_PntFlowUniform,max_LvDiff,max_PntOvrlp
INTEGER lv_Uniform,ratio_FlowWall
INTEGER flag_Add1Lyr
REAL(8) l_bc,y_goal,y1,ER1,wall_Position,y_goal_input
REAL(8) tmp_Real,tmp_Real1,tmp_Real2
INTEGER,DIMENSION(:),ALLOCATABLE ::  n_GrdLyr,n_PntFlow,lv_Grid
REAL(8),DIMENSION(:),ALLOCATABLE ::  y1_Lyr,ER_Lyr,yn_Lyr,x1_Lyr
REAL(8),DIMENSION(:),ALLOCATABLE ::  thcknss_Lyr,thcknss_TotDue
CHARACTER(LEN=20) FILENAME
20 FORMAT(A50)
21 FORMAT(A50,1X,I8)
22 FORMAT(A50,1X,I8,1X,I8)
31 FORMAT(A50,1X,F10.8)
!--------------------------------------------------------------------------!
WRITE(FILENAME,'(A15)') "file_WallBC_Gen"
OPEN(10,FILE=FILENAME,STATUS="OLD")
READ(10,*) 
READ(10,*) y1
READ(10,*) 
READ(10,*) ER1
READ(10,*) 
READ(10,*) ratio_FlowWall
READ(10,*) 
READ(10,*) l_bc
READ(10,*) 
READ(10,*) y_goal
READ(10,*) 
READ(10,*) lv_Uniform
READ(10,*) 
READ(10,*) n_PntOvrlp
READ(10,*) 
READ(10,*) wall_Position
CLOSE(10)
y_goal_input=y_goal
WRITE(*,31) "y1=:",y1
WRITE(*,31) "ER1=:",ER1
WRITE(*,21) "ratio_FlowWall=:",ratio_FlowWall
WRITE(*,31) "l_bc=:",l_bc
WRITE(*,31) "y_goal_input=:",y_goal_input
WRITE(*,21) "lv_Uniform=:",lv_Uniform
WRITE(*,21) "n_PntOvrlp=:",n_PntOvrlp
WRITE(*,31) "wall_Position=:",wall_Position
n_PntFlowUniform=l_bc/y_goal
CALL nice_num(n_PntFlowUniform,lv_Uniform)
WRITE(*,'(A60)') "------------------------------------------------------&
                 --------------------"
WRITE(*,21) "n_PntFlowUniform(Origin)=:",n_PntFlowUniform
y_goal=l_bc/n_PntFlowUniform
WRITE(*,31) "y_goal=:",y_goal
!--------------------------------------------------------------------------!
n_LvDown=0
tmp_Int=0
flag_Add1Lyr=0
DO WHILE(tmp_Int<=ratio_FlowWall)
  n_LvDown=n_LvDown+1
  tmp_Int=2**n_LvDown
ENDDO
n_LvDown=n_LvDown-1
max_LvDiff=n_LvDown
max_PntOvrlp=n_PntOvrlp*2**max_LvDiff
WRITE(*,21) "max_PntOvrlp=:",max_PntOvrlp
WRITE(*,21) "n_LvDown/max_LvDiff=:",n_LvDown
!--------------------------------------------------------------------------!
ratio_ynVSy1=y_goal/y1
WRITE(*,21) "ratio_ynVSy1=:",ratio_ynVSy1
n_LvTot=0
tmp_Int=0
DO WHILE(tmp_Int<=ratio_ynVSy1)
  n_LvTot=n_LvTot+1
  tmp_Int=2**n_LvTot
ENDDO
n_LvTot=n_LvTot-1
WRITE(*,21) "n_LvTot=:",n_LvTot
n_BlkLyr=n_LvTot/n_LvDown
200 IF(flag_Add1Lyr==1) n_BlkLyr=n_BlkLyr+1
WRITE(*,21) "n_BlkLyr=:",n_BlkLyr
!--------------------------------------------------------------------------!
ALLOCATE(n_GrdLyr(n_BlkLyr),n_PntFlow(n_BlkLyr),lv_Grid(n_BlkLyr))
ALLOCATE(y1_Lyr(n_BlkLyr),ER_Lyr(n_BlkLyr))
ALLOCATE(yn_Lyr(n_BlkLyr),x1_Lyr(n_BlkLyr))
ALLOCATE(thcknss_Lyr(n_BlkLyr),thcknss_TotDue(n_BlkLyr))
thcknss_TotDue=0.d0
DO I=1,n_BlkLyr
  WRITE(*,'(A60)') "------------------------------------------------------&
                   --------------------"
  WRITE(*,21) "i_BlkLyr=:",I
  IF(I==1) THEN
    y1_Lyr(I)=y1
    x1_Lyr(I)=y1*ratio_FlowWall
    n_PntFlow(I)=l_bc/x1_Lyr(I)
    CALL nice_num(n_PntFlow(1),n_LvTot)
    ER_Lyr(I)=ER1
    lv_Grid(I)=n_LvTot-n_LvDown
  ELSEIF(I<n_BlkLyr) THEN
    y1_Lyr(I)=yn_Lyr(I-1)*ER_Lyr(I-1)
    n_PntFlow(I)=(n_PntFlow(I-1)-1)/2**n_LvDown
    ER_Lyr(I)=ER_Lyr(I-1)
    lv_Grid(I)=lv_Grid(I-1)-n_LvDown
  ELSE
    IF(flag_Add1Lyr==1) n_LvDown=n_LvDown-1
    y1_Lyr(I)=yn_Lyr(I-1)*ER_Lyr(I-1)
    n_PntFlow(I)=(n_PntFlow(I-1)-1)/2**n_LvDown
    ER_Lyr(I)=ER_Lyr(I-1)
    lv_Grid(I)=lv_Grid(I-1)-n_LvDown
  ENDIF
  WRITE(*,21) "lv_Grid:",lv_Grid(I)
  WRITE(*,31) "y1_Lyr:",y1_Lyr(I)
  x1_Lyr(I)=l_bc/n_PntFlow(I)
  yn_Lyr(I)=x1_Lyr(I)
  WRITE(*,31) "x1_Lyr/yn_Lyr=:",x1_Lyr(I)
  WRITE(*,22) "n_PntFlow:",n_PntFlow(I)
  n_GrdLyr(I)=DLOG(yn_Lyr(I)/y1_Lyr(I))/DLOG(ER_Lyr(I))+1
!  WRITE(*,21) "n_GrdLyr:",n_GrdLyr(I)
  CALL nice_num(n_GrdLyr(I),lv_Uniform)
  WRITE(*,21) "n_GrdLyr:",n_GrdLyr(I)
  IF(n_GrdLyr(I)<max_PntOvrlp) THEN
    WRITE(*,21) "Not enough grid layers in wall block:",I
    STOP 1
  ENDIF
  tmp_Real1=yn_Lyr(I)/y1_Lyr(I)
  tmp_Real2=1.d0/(REAL(n_GrdLyr(I))-1.d0)
  ER_Lyr(I)=tmp_Real1**tmp_Real2
  WRITE(*,31) "ER_Lyr:",ER_Lyr(I)
  thcknss_Lyr(I)=y1_Lyr(I)*(1-ER_Lyr(I)**n_GrdLyr(I))/(1-ER_Lyr(I))
!  WRITE(*,31) "thcknss_Lyr:",thcknss_Lyr(I)
  IF(I==1) THEN
    thcknss_TotDue(I)=thcknss_Lyr(I)
  ELSE
    thcknss_TotDue(I)=thcknss_TotDue(I-1)+thcknss_Lyr(I)
  ENDIF
  WRITE(*,31) "thcknss_TotDue:",thcknss_TotDue(I)
ENDDO
IF(y_goal/x1_Lyr(n_BlkLyr)>2 .AND. n_LvDown>1) THEN
  WRITE(*,'(A60)') "------------------------------------------------------&
                   --------------------"
  WRITE(*,'(A60)') "######################################################&
                   ####################"
  WRITE(*,20) "Set one more block layer."
  flag_Add1Lyr=1
  DEALLOCATE(n_GrdLyr,n_PntFlow,y1_Lyr,ER_Lyr)
  DEALLOCATE(yn_Lyr,x1_Lyr,thcknss_Lyr,lv_Grid,thcknss_TotDue)
  GOTO 200
ENDIF
tmp_Real=y_goal
y_goal=yn_Lyr(n_BlkLyr)
n_PntFlowUniform=l_bc/y_goal
CALL nice_num(n_PntFlowUniform,lv_Uniform)
WRITE(FILENAME,'(A18)') "file_WallBC_Output"
OPEN(20,FILE=FILENAME,STATUS="UNKNOWN")
WRITE(20,31) "y1=:",y1
WRITE(20,31) "ER1=:",ER1
WRITE(20,21) "ratio_FlowWall=:",ratio_FlowWall
WRITE(20,31) "l_bc=:",l_bc
WRITE(20,31) " y_goal_input=:",y_goal_input
WRITE(20,31) "y_goal_origin=:",tmp_Real
WRITE(20,31) "   y_goal_now=:",y_goal
WRITE(20,21) "lv_Uniform=:",lv_Uniform
WRITE(20,21) "n_PntOvrlp=:",n_PntOvrlp
WRITE(20,21) "n_PntFlowUniform=:",n_PntFlowUniform
WRITE(20,31) "wall_Position=:",wall_Position
WRITE(20,'(A65)') "------------------------------------------------------&
                 --------------------"
DO I=1,n_BlkLyr
  WRITE(20,21) "Index of  block,I=:",I
  WRITE(20,21) "Grid level in this block,lv_Grid=:",(-1)*lv_Grid(I)
  WRITE(20,31) "First layer thickness in this block,y1=:",y1_Lyr(I)
  WRITE(20,31) "Last layer thickness in this block,yn=:",yn_Lyr(I)
  WRITE(20,21) "No. points in flow direction,n_PntFlow=:",n_PntFlow(I)
  WRITE(20,21) "No. points in wall normal direction,n_GrdLyr=:",n_GrdLyr(I)
  WRITE(20,31) "Expansion ratio in this block,ER=:",ER_Lyr(I)
  WRITE(20,31) "Total thickness due in this block,thcknss_TotDue=:",thcknss_TotDue(I)
  WRITE(20,31) "Block boundary position,wall-=:",wall_Position-thcknss_TotDue(I)
  WRITE(20,31) "Block boundary position,wall+=:",wall_Position+thcknss_TotDue(I)
  WRITE(20,'(A65)') "------------------------------------------------------&
                   --------------------"
ENDDO
CLOSE(20)

END

SUBROUTINE nice_num(n,lv)
IMPLICIT NONE
INTEGER n,lv
n=((n-1)/2**lv+1)*2**lv+1
RETURN
END SUBROUTINE


