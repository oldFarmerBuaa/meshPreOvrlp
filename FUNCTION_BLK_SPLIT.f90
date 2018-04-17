!-----------------------------------------------------------------------!
!   Function of block splitting                                         !
!   Split previous block into two sub-blocks				!
!   1) i_BlkSplit = index of block to be split				!
!   2) i_Split    = index of split out of total split in a block	!
!   3) i_BlkGlb   = global block index by now before final split	!
!   4) J,K        = block index and super-block index			!
!-----------------------------------------------------------------------!
subroutine function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
!
USE READ_GRIDCON_Pre
USE system_var
USE block_overlaping
USE block_splitting
!
IMPLICIT NONE
! Block splitting
INTEGER I,J,K,i_I,i_J,i_K
INTEGER i_BlkSplit,i_Split,i_BlkGlb,i_BlkSplit1,i_BlkSplit2
INTEGER KK_PreStp,JJ_PreStp,II_PreStp,n_PntPreStp,max_IJK_Loc
INTEGER IJ_start,IJ_end
INTEGER n_PntNowBlk1,n_PntNowBlk2
INTEGER st_I1,st_J1,st_K1,ed_I1,ed_J1,ed_K1
INTEGER st_I2,st_J2,st_K2,ed_I2,ed_J2,ed_K2
INTEGER mid_II,mid_JJ,mid_KK,n_PntHlfOvrlp
INTEGER II_NowBlk_1,JJ_NowBlk_1,KK_NowBlk_1,II_NowBlk_2,JJ_NowBlk_2,KK_NowBlk_2
INTEGER i_PntPre,i_PntNow1,i_PntNow2
INTEGER II,JJ,KK
INTEGER isidtp_TmpPre(6),isidtp_TmpAft(6,2),isidtp_Aft(6)
INTEGER Lv_OffUnfrmTmpPre(3),Lv_OffUnfrmTmpAft(3,2),Lv_OffUnfrmAft(3)
REAL(8),DIMENSION(:),ALLOCATABLE :: X_PreStp,Y_PreStp,Z_PreStp
REAL(8),DIMENSION(:),ALLOCATABLE :: X_NowBlk_1,Y_NowBlk_1,Z_NowBlk_1
REAL(8),DIMENSION(:),ALLOCATABLE :: X_NowBlk_2,Y_NowBlk_2,Z_NowBlk_2

20 FORMAT(A50)
21 FORMAT(A50,1X,I8)
22 FORMAT(A50,1X,I8,1X,I8)
23 FORMAT(A50,1X,I8,1X,I8,1X,I8)
!-----------------------------------------------------------------------!
! Reading interval block to be split
IF(i_Split==1) THEN
  II_PreStp=II_toSplit(J,K)
  JJ_PreStp=JJ_toSplit(J,K)
  KK_PreStp=KK_toSplit(J,K)
  n_PntPreStp=II_PreStp*JJ_PreStp*KK_PreStp
  ALLOCATE(X_PreStp(n_PntPreStp))
  ALLOCATE(Y_PreStp(n_PntPreStp))
  ALLOCATE(Z_PreStp(n_PntPreStp))
  X_PreStp(1:n_PntPreStp)=X_toSplit(IJ_start:IJ_end,K)
  Y_PreStp(1:n_PntPreStp)=Y_toSplit(IJ_start:IJ_end,K)
  Z_PreStp(1:n_PntPreStp)=Z_toSplit(IJ_start:IJ_end,K)
  isidtp_TmpPre(1:6)=isidtp_Pre(1:6,J,K)
  Lv_OffUnfrmTmpPre(1:3)=Lv_OffUniform(1:3,J,K)
!  WRITE(*,*)"X_PreStp(1):",X_PreStp(1)
!  WRITE(*,*)"X_PreStp(end):",X_PreStp(n_PntPreStp)
ELSE
  WRITE(FILENAME,'(A13,I8.8)')"GRIDBLO_INTVL",i_BlkSplit
  KK_PreStp=1
  OPEN(30,FILE=FILENAME,STATUS="OLD")
  IF(cas_Dim=="2D") THEN
    READ(30,*)
    READ(30,*) (Lv_OffUnfrmTmpPre(I),I=1,3)
    READ(30,*)
    READ(30,*) (isidtp_TmpPre(I),I=1,6)
    READ(30,*)
    READ(30,*) JJ_PreStp,II_PreStp
    n_PntPreStp=II_PreStp*JJ_PreStp*KK_PreStp
    ALLOCATE(X_PreStp(n_PntPreStp))
    ALLOCATE(Y_PreStp(n_PntPreStp))
    ALLOCATE(Z_PreStp(n_PntPreStp))
    READ(30,*) (X_PreStp(i_I),i_I=1,n_PntPreStp)
    READ(30,*) (Y_PreStp(i_I),i_I=1,n_PntPreStp)
    READ(30,*) (Z_PreStp(i_I),i_I=1,n_PntPreStp)
!    WRITE(*,*)"X_PreStp(1):",X_PreStp(1)
!    WRITE(*,*)"X_PreStp(end):",X_PreStp(n_PntPreStp)
  ELSE
    READ(30,*)
    READ(30,*) (Lv_OffUnfrmTmpPre(I),I=1,3)
    READ(30,*)
    READ(30,*) (isidtp_TmpPre(I),I=1,6)
    READ(30,*)
    READ(30,*) KK_PreStp,JJ_PreStp,II_PreStp
    n_PntPreStp=II_PreStp*JJ_PreStp*KK_PreStp
    ALLOCATE(X_PreStp(n_PntPreStp))
    ALLOCATE(Y_PreStp(n_PntPreStp))
    ALLOCATE(Z_PreStp(n_PntPreStp))
    READ(30,*) (X_PreStp(i_I),i_I=1,n_PntPreStp)
    READ(30,*) (Y_PreStp(i_I),i_I=1,n_PntPreStp)
    READ(30,*) (Z_PreStp(i_I),i_I=1,n_PntPreStp)
  ENDIF
  CLOSE(30)
ENDIF
!-----------------------------------------------------------------------!
! Split previous block into two sub-blocks
max_IJK_Loc=MAX(II_PreStp,JJ_PreStp,KK_PreStp)
n_PntHlfOvrlp=n_PntOvrlpInSplit/2
!Assign start and end II/JJ/KK value of each sub block
!Firstly, assign II/JJ/KK to fully cover PreStp block
st_I1=1
st_J1=1
st_K1=1
ed_I1=II_PreStp
ed_J1=JJ_PreStp
ed_K1=KK_PreStp
st_I2=1
st_J2=1
st_K2=1
ed_I2=II_PreStp
ed_J2=JJ_PreStp
ed_K2=KK_PreStp
!Secondly, adjust value by cutting position, that is , along II/JJ or KK.    
!Cut on the edge with the most points.
isidtp_TmpAft(1:6,1)=isidtp_TmpPre(1:6)
isidtp_TmpAft(1:6,2)=isidtp_TmpPre(1:6)
Lv_OffUnfrmTmpAft(1:3,1)=Lv_OffUnfrmTmpPre(1:3)
Lv_OffUnfrmTmpAft(1:3,2)=Lv_OffUnfrmTmpPre(1:3)
IF(max_IJK_Loc==II_PreStp) THEN
!  WRITE(*,20) "Split by I."
  mid_II=II_PreStp/2+1
  II_NowBlk_1=mid_II+n_PntHlfOvrlp
  JJ_NowBlk_1=JJ_PreStp
  KK_NowBlk_1=KK_PreStp
  II_NowBlk_2=II_PreStp-mid_II+1+n_PntHlfOvrlp
  JJ_NowBlk_2=JJ_PreStp
  KK_NowBlk_2=KK_PreStp
  ed_I1=II_NowBlk_1
  st_I2=mid_II-n_PntHlfOvrlp
  isidtp_TmpAft(2,1)=9000
  isidtp_TmpAft(1,2)=9000
ELSEIF(max_IJK_Loc==JJ_PreStp) THEN
!  WRITE(*,20) "Split by J."
  mid_JJ=JJ_PreStp/2+1
  II_NowBlk_1=II_PreStp
  JJ_NowBlk_1=mid_JJ+n_PntHlfOvrlp
  KK_NowBlk_1=KK_PreStp
  II_NowBlk_2=II_PreStp
  JJ_NowBlk_2=JJ_PreStp-mid_JJ+1+n_PntHlfOvrlp
  KK_NowBlk_2=KK_PreStp
  ed_J1=JJ_NowBlk_1
  st_J2=mid_JJ-n_PntHlfOvrlp
  isidtp_TmpAft(4,1)=9000
  isidtp_TmpAft(3,2)=9000
ELSE
!  WRITE(*,20) "Split by K."
  mid_KK=KK_PreStp/2+1
  II_NowBlk_1=II_PreStp
  JJ_NowBlk_1=JJ_PreStp
  KK_NowBlk_1=mid_KK+n_PntHlfOvrlp
  II_NowBlk_2=II_PreStp
  JJ_NowBlk_2=JJ_PreStp
  KK_NowBlk_2=KK_PreStp-mid_KK+1+n_PntHlfOvrlp
  ed_K1=JJ_NowBlk_1
  st_K2=mid_KK-n_PntHlfOvrlp
  isidtp_TmpAft(6,1)=9000
  isidtp_TmpAft(5,2)=9000
ENDIF
n_PntNowBlk1=II_NowBlk_1*JJ_NowBlk_1*KK_NowBlk_1
n_PntNowBlk2=II_NowBlk_2*JJ_NowBlk_2*KK_NowBlk_2
!WRITE(*,22) "st_I1,ed_I1:",st_I1,ed_I1
!WRITE(*,22) "st_I2,ed_I2:",st_I2,ed_I2
!WRITE(*,22) "st_J1,ed_J1:",st_J1,ed_J1
!WRITE(*,22) "st_J2,ed_J2:",st_J2,ed_J2
!WRITE(*,23) "I/J/K_NowBlk1:",II_NowBlk_1,JJ_NowBlk_1,KK_NowBlk_1
!WRITE(*,23) "I/J/K_NowBlk2:",II_NowBlk_2,JJ_NowBlk_2,KK_NowBlk_2
!WRITE(*,22) "n_PntNowBlk_1/2:",n_PntNowBlk1,n_PntNowBlk2
ALLOCATE(X_NowBlk_1(n_PntNowBlk1))
ALLOCATE(Y_NowBlk_1(n_PntNowBlk1))
ALLOCATE(Z_NowBlk_1(n_PntNowBlk1))
ALLOCATE(X_NowBlk_2(n_PntNowBlk2))
ALLOCATE(Y_NowBlk_2(n_PntNowBlk2))
ALLOCATE(Z_NowBlk_2(n_PntNowBlk2))
DO i_K=st_K1,ed_K1
  DO i_J=st_J1,ed_J1
    DO i_I=st_I1,ed_I1
      i_PntPre=(i_K-1)*JJ_PreStp*II_PreStp+(i_J-1)*II_PreStp+i_I
      i_PntNow1=(i_K-1)*JJ_NowBlk_1*II_NowBlk_1+(i_J-1)*II_NowBlk_1+i_I
!      WRITE(*,*) "i_PntPre,i_PntNow1:",i_PntPre,i_PntNow1
      X_NowBlk_1(i_PntNow1)=X_PreStp(i_PntPre)
      Y_NowBlk_1(i_PntNow1)=Y_PreStp(i_PntPre)
      Z_NowBlk_1(i_PntNow1)=Z_PreStp(i_PntPre)
    ENDDO
  ENDDO
ENDDO    
DO i_K=st_K2,ed_K2
  DO i_J=st_J2,ed_J2
    DO i_I=st_I2,ed_I2
      i_PntPre=(i_K-1)*JJ_PreStp*II_PreStp+(i_J-1)*II_PreStp+i_I
      IF(max_IJK_Loc==II_PreStp) THEN
      i_PntNow2=(i_K-1)*JJ_NowBlk_2*II_NowBlk_2+(i_J-1)*II_NowBlk_2+i_I-st_I2+1
      ELSEIF(max_IJK_Loc==JJ_PreStp) THEN
      i_PntNow2=(i_K-1)*JJ_NowBlk_2*II_NowBlk_2+(i_J-st_J2+1-1)*II_NowBlk_2+i_I
      ELSE
      i_PntNow2=(i_K-st_K2+1-1)*JJ_NowBlk_2*II_NowBlk_2+(i_J-1)*II_NowBlk_2+i_I
      ENDIF
      X_NowBlk_2(i_PntNow2)=X_PreStp(i_PntPre)
      Y_NowBlk_2(i_PntNow2)=Y_PreStp(i_PntPre)
      Z_NowBlk_2(i_PntNow2)=Z_PreStp(i_PntPre)            
    ENDDO
  ENDDO
ENDDO
!-----------------------------------------------------------------------!
! Output interval and finale GRIDBLOCKS
IF(i_Split==n_Split(J,K)) THEN
! Last split, output final gird block
  II=II_NowBlk_1
  JJ=JJ_NowBlk_1
  KK=KK_NowBlk_1
  i_BlkGlb=i_BlkGlb+1
  isidtp_Aft(1:6)=isidtp_TmpAft(1:6,1)
  Lv_OffUnfrmAft(1:3)=Lv_OffUnfrmTmpAft(1:3,1)
  CALL  function_Output_Grid(II,JJ,KK,i_BlkGlb,&
        X_NowBlk_1,Y_NowBlk_1,Z_NowBlk_1,&
        n_PntNowBlk1,isidtp_Aft,Lv_OffUnfrmAft)
  II=II_NowBlk_2
  JJ=JJ_NowBlk_2
  KK=KK_NowBlk_2
  i_BlkGlb=i_BlkGlb+1
  isidtp_Aft(1:6)=isidtp_TmpAft(1:6,2)
  Lv_OffUnfrmAft(1:3)=Lv_OffUnfrmTmpAft(1:3,2)
  CALL  function_Output_Grid(II,JJ,KK,i_BlkGlb,&
        X_NowBlk_2,Y_NowBlk_2,Z_NowBlk_2,&
        n_PntNowBlk2,isidtp_Aft,Lv_OffUnfrmAft)
ELSE
  i_BlkSplit1=i_BlkSplit*10+1
  i_BlkSplit2=i_BlkSplit*10+2
! Output sub-block 1a
  100 FORMAT(3F16.12)
  WRITE(FILENAME,'(A13,I8.8)')"GRIDBLO_INTVL",i_BlkSplit1
  OPEN(30,FILE=FILENAME,STATUS="UNKNOWN")
  IF(cas_Dim=="2D") THEN
    WRITE(*,*) "Writing file1:",FILENAME
    WRITE(30,*) "Grid level is :Lv_OffUniform(1:3):"
    WRITE(30,*) (Lv_OffUnfrmTmpAft(I,1),I=1,3)
    WRITE(30,*) "Boundary Conditions:isidtp(1:6):"
    WRITE(30,*) (isidtp_TmpAft(I,1),I=1,6)
    WRITE(30,*) "JJ,II and X(1:n),Y(1:n),Z(1:n):"
    WRITE(30,*) JJ_NowBlk_1,II_NowBlk_1
    WRITE(30,100) (X_NowBlk_1(i_I),i_I=1,n_PntNowBlk1)
    WRITE(30,100) (Y_NowBlk_1(i_I),i_I=1,n_PntNowBlk1)
    WRITE(30,100) (Z_NowBlk_1(i_I),i_I=1,n_PntNowBlk1)
  ELSE
    WRITE(*,*) "Writing file1:",FILENAME
    WRITE(30,*) "Grid level is :Lv_OffUniform(1:3):"
    WRITE(30,*) (Lv_OffUnfrmTmpAft(I,1),I=1,3)
    WRITE(30,*) "Boundary Conditions:isidtp(1:6):"
    WRITE(30,*) (isidtp_TmpAft(I,1),I=1,6)
    WRITE(30,*) "KK,JJ,II and X(1:n),Y(1:n),Z(1:n):"
    WRITE(30,*) KK_NowBlk_1,JJ_NowBlk_1,II_NowBlk_1
    WRITE(30,100) (X_NowBlk_1(i_I),i_I=1,n_PntNowBlk1)
    WRITE(30,100) (Y_NowBlk_1(i_I),i_I=1,n_PntNowBlk1)
    WRITE(30,100) (Z_NowBlk_1(i_I),i_I=1,n_PntNowBlk1)
  ENDIF
  CLOSE(30)
!Output sub-block 2
  WRITE(FILENAME,'(A13,I8.8)')"GRIDBLO_INTVL",i_BlkSplit2
  OPEN(30,FILE=FILENAME,STATUS="UNKNOWN")
  IF(cas_Dim=="2D") THEN
    WRITE(*,*) "Writing file2:",FILENAME
    WRITE(30,*) "Grid level is :Lv_OffUniform(1:3):"
    WRITE(30,*) (Lv_OffUnfrmTmpAft(I,2),I=1,3)
    WRITE(30,*) "Boundary Conditions:isidtp(1:6):"
    WRITE(30,*) (isidtp_TmpAft(I,2),I=1,6)
    WRITE(30,*) "JJ,II and X(1:n),Y(1:n),Z(1:n):"
    WRITE(30,*) JJ_NowBlk_2,II_NowBlk_2
    WRITE(30,100) (X_NowBlk_2(i_I),i_I=1,n_PntNowBlk2)
    WRITE(30,100) (Y_NowBlk_2(i_I),i_I=1,n_PntNowBlk2)
    WRITE(30,100) (Z_NowBlk_2(i_I),i_I=1,n_PntNowBlk2)
  ELSE
    WRITE(*,*) "Writing file2:",FILENAME
    WRITE(30,*) "Grid level is :Lv_OffUniform(1:3):"
    WRITE(30,*) (Lv_OffUnfrmTmpAft(I,2),I=1,3)
    WRITE(30,*) "Boundary Conditions:isidtp(1:6):"
    WRITE(30,*) (isidtp_TmpAft(I,2),I=1,6)
    WRITE(30,*) "KK,JJ,II and X(1:n),Y(1:n),Z(1:n):"
    WRITE(30,*) KK_NowBlk_2,JJ_NowBlk_2,II_NowBlk_2
    WRITE(30,100) (X_NowBlk_2(i_I),i_I=1,n_PntNowBlk2)
    WRITE(30,100) (Y_NowBlk_2(i_I),i_I=1,n_PntNowBlk2)
    WRITE(30,100) (Z_NowBlk_2(i_I),i_I=1,n_PntNowBlk2)
  ENDIF
  CLOSE(30)
ENDIF
RETURN
END

!-----------------------------------------------------------------------!
!   Output GRID blocks in splitting                                     !
!-----------------------------------------------------------------------!
SUBROUTINE function_Output_Grid(II,JJ,KK,i_BlkGlb,X,Y,Z,n_Point,isidtp_Output,Lv_OffUnfrmOutput)
!
USE READ_GRIDCON_Pre
!
IMPLICIT NONE
INTEGER I,J,K,i_IJK,i_JIK
INTEGER II,JJ,KK,i_BlkGlb,n_Point
INTEGER isidtp_Output(6),Lv_OffUnfrmOutput(3)
REAL(8) X(n_Point),Y(n_Point),Z(n_Point)
REAL(8) X_JIK(n_Point),Y_JIK(n_Point),Z_JIK(n_Point)
CHARACTER(LEN=20) FILENAME1
CHARACTER(LEN=50) FILENAME2

21 FORMAT(A50,1X,I8)
23 FORMAT(A50,1X,I8,1X,I8,1X,I8)
! Second final output
WRITE(FILENAME1,'(A10,I4.4)')"GRIDBLO_BC",i_BlkGlb
! Final output in which grid point is arranged by 
! i_Point=(K-1)*IMAX*JMAX+(I-1)*JMAX+J, TRANSFORM I-J-K TO J-I-K COUNT
! 2DDNS CODE IS J-I-K
WRITE(FILENAME2,'(A19,I3.3)')"output_GRID/GRIDBLO",i_BlkGlb
IF(i_BlkGLb>999) CALL PSEXIT("n_BlkGlb>999, increase block name digit.")
OPEN(30,FILE=FILENAME1,STATUS="UNKNOWN")
IF(cas_Dim=="2D") THEN
  IF(II<17 .OR. JJ<17) THEN
    WRITE(*,21) "Warning! Less than 17 points in block:",i_BlkGlb
  ELSEIF(II<7 .OR. JJ<7) THEN
    WRITE(*,21) "Idx of block:",i_BlkGlb
    CALL PSEXIT("Error! Less than 7 points in some blocks.")
  ENDIF
  WRITE(30,*) "Grid level is :Lv_OffUniform(1:3):"
  WRITE(30,*) (Lv_OffUnfrmOutput(I),I=1,3)
  WRITE(30,*) "Boundary Conditions:isidtp(1:6):"
  WRITE(30,*) (isidtp_Output(I),I=1,6)
  WRITE(30,*) "JJ,II and X(1:n),Y(1:n):"
  WRITE(30,*) JJ,II
  WRITE(*,23) "Index of block after splitting:",i_BlkGlb
  WRITE(*,23) "Total points in this block:",n_Point
!WRITE(*,*) "--------------------------------------&
!            ----------------------------------------"
  WRITE(30,100) (X(I),I=1,n_Point)
  WRITE(30,100) (Y(I),I=1,n_Point)
ELSE
  IF(II<17 .OR. JJ<17 .OR. KK<17) THEN
    WRITE(*,21) "Warning! Less than 17 points in block:",i_BlkGlb
  ELSEIF(II<7 .OR. JJ<7 .OR. KK<7) THEN
    WRITE(*,21) "Idx of block:",i_BlkGlb
    CALL PSEXIT("Error! Less than 7 points in some blocks.")
  ENDIF
  WRITE(30,*) "Grid level is :Lv_OffUniform(1:3):"
  WRITE(30,*) (Lv_OffUnfrmOutput(I),I=1,3)
  WRITE(30,*) "Boundary Conditions:isidtp(1:6):"
  WRITE(30,*) (isidtp_Output(I),I=1,6)
  WRITE(30,*) "KK,JJ,II and X(1:n),Y(1:n),Z(1:n):"
  WRITE(30,*) KK,JJ,II
  WRITE(30,100) (X(I),I=1,n_Point)
  WRITE(30,100) (Y(I),I=1,n_Point)
  WRITE(30,100) (Z(I),I=1,n_Point)
ENDIF
CLOSE(30)
DO K=1,KK
  DO J=1,JJ
    DO I=1,II
      i_JIK=(K-1)*II*JJ+(I-1)*JJ+J
      i_IJK=(K-1)*II*JJ+(J-1)*II+I
      X_JIK(i_JIK)=X(i_IJK)
      Y_JIK(i_JIK)=Y(i_IJK)
      Z_JIK(i_JIK)=Z(i_IJK)
    ENDDO
  ENDDO
ENDDO
OPEN(30,FILE=FILENAME2,STATUS="UNKNOWN")
IF(cas_Dim=="2D") THEN
  WRITE(30,*) II,JJ
  WRITE(30,100) (X_JIK(i_JIK),i_JIK=1,n_Point)
  WRITE(30,100) (Y_JIK(i_JIK),i_JIK=1,n_Point)
ELSE
  WRITE(30,*) II,JJ,KK
  WRITE(30,100) (X_JIK(i_JIK),i_JIK=1,n_Point)
  WRITE(30,100) (Y_JIK(i_JIK),i_JIK=1,n_Point)
  WRITE(30,100) (Z_JIK(i_JIK),i_JIK=1,n_Point)
ENDIF
CLOSE(30)
100 FORMAT(3F16.12)
RETURN
END SUBROUTINE
       
