!-----------------------------------------------------------------------!
!    Read grid from mesh file						!
!-----------------------------------------------------------------------!
subroutine read_GrdAndAllocate
!
USE read_gridcon_pre
USE read_mesh_original
USE system_var
USE point_leaping
USE block_overlaping
USE block_splitting
!
IMPLICIT NONE
!-----------------------------------------------------------------------!
! Point index
INTEGER IJ_start,IJ_end
INTEGER i_I,i_J,i_K,I,J,K,III,JJJ,KKK,i_Point0,i_PointLP,&
        i_Point,i_PointLoca,i_Point2,i_PointLoc,i_PointLocPre,&
        i_BlkGlb,i_BlkGlb2
INTEGER I_Leap,J_Leap,K_Leap
INTEGER i_Loc,i_I_aft,i_J_aft,i_K_aft,i_Loop,j_Loop,k_Loop,i_l,n_PointAdd
INTEGER step_PntLp(3)
REAL(8),DIMENSION(:),ALLOCATABLE   :: X_Loc,Y_Loc,Z_Loc
REAL(8),DIMENSION(:),ALLOCATABLE   :: X_LocPre,Y_LocPre,Z_LocPre

! Format of screen/logfile output, 2N = N data output
20 FORMAT(A50)
21 FORMAT(A50,1X,I8)
22 FORMAT(A50,1X,I8,1X,I8)
23 FORMAT(A50,1X,I8,1X,I8,1X,I8)
100   FORMAT(2F16.12)
!-----------------------------------------------------------------------!
! Start
WRITE(LOGFILE,'(A18)') "meshProcessing.log"
OPEN(19,FILE=LOGFILE,STATUS="UNKNOWN")
OPEN(10,FILE="GRIDCON_Pre",STATUS="OLD")
!-----------------------------------------------------------------------!
!   Read GRIDCON_Pre for total blocks in each sup-block			!
!-----------------------------------------------------------------------!
READ(10,*) 
READ(10,*) cas_Dim
READ(10,*) 
READ(10,*) mesh_Type
READ(10,*) 
READ(10,*) n_BlkSupPre 
WRITE(*,21) "Number of super blocks: ",n_BlkSupPre
READ(10,*) 
READ(10,*) 
!-----------------------------------------------------------------------!
! Allocate global arrays
ALLOCATE(                      n_BlkPre(n_BlkSupPre))
n_BlkGlbPre=0
DO K=1,n_BlkSupPre
  READ(10,*) tmp_Int,n_BlkPre(K)
  n_BlkGlbPre=n_BlkGlbPre+n_BlkPre(K)
ENDDO
max_BlkGlbPre=MAXVAL(n_BlkPre)
! BC values,1:6 for 3D
! Lv_OffUniform = Lv_OffSet + Lv_Coarsen
ALLOCATE(                    n_PointPre(n_BlkSupPre))
ALLOCATE(                    n_PntAftLp(n_BlkSupPre))
ALLOCATE(                  n_PntToOvrlp(n_BlkSupPre))
ALLOCATE(n_PntAftOvrlpBfBlk(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(      n_PntInBlk(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(      n_PntBfBlk(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(    n_toSplitLoc(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(                 n_PntAftOvrlp(n_BlkSupPre))
ALLOCATE(   n_PointLocPre(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(       n_PointLP(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(   n_PntOvrlpLoc(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(              max_GlueInBlkPre(n_BlkSupPre))
ALLOCATE(                        II_Blk(n_BlkSupPre))
ALLOCATE(                        JJ_Blk(n_BlkSupPre))
ALLOCATE(                        KK_Blk(n_BlkSupPre))
ALLOCATE(    isidtp_Pre(6,max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(     dir_Swipe(6,max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(                  Lv_Uniform(3,n_BlkSupPre))
ALLOCATE(     Lv_OffSet(3,max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(    Lv_Coarsen(3,max_BlkGlbPre,n_BlkSupPre))
ALLOCATE( Lv_OffUniform(3,max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(Lv_OffUniformTot(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(  Lv_DiffVsMin(3,max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(          II_Pre(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(          JJ_Pre(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(          KK_Pre(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(         II_toLP(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(         JJ_toLP(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(         KK_toLP(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(        II_PntLp(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(        JJ_PntLp(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(        KK_PntLp(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(      II_toOvrlp(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(      JJ_toOvrlp(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(      KK_toOvrlp(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(     II_PntOvrlp(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(     JJ_PntOvrlp(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(     KK_PntOvrlp(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(      II_toSplit(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(      JJ_toSplit(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(      KK_toSplit(max_BlkGlbPre,n_BlkSupPre))
!-----------------------------------------------------------------------!
!   Read GRIDCON_Pre for BCs in each block				!
!-----------------------------------------------------------------------!
READ(10,*)
DO I=1,n_BlkSupPre
  READ(10,*)
  READ(10,*)
  READ(10,*) Lv_Uniform(1,I),Lv_Uniform(2,I),Lv_Uniform(3,I)
ENDDO
!-----------------------------------------------------------------------!
!#14- No. Level of global coarsen, Lv_CoarsenGlb
READ(10,*)
READ(10,*) Lv_CoarsenGlb
DO I=1,11
  READ(10,*)
ENDDO
DO K=1,n_BlkSupPre
  READ(10,*) 
  DO J=1,n_BlkPre(K)
    READ(10,*) tmp_Int,isidtp_Pre(1,J,k),isidtp_Pre(2,J,k),&
               isidtp_Pre(3,J,k),isidtp_Pre(4,J,k),&
               isidtp_Pre(5,J,k),isidtp_Pre(6,J,k),&
               Lv_OffSet(1,J,K),Lv_OffSet(2,J,K),Lv_OffSet(3,J,K),&
               Lv_Coarsen(1,J,K),Lv_Coarsen(2,J,K),Lv_Coarsen(3,J,K)
    DO I=1,3
    Lv_OffUniform(I,J,K) = Lv_OffSet(I,J,K) + Lv_Coarsen(I,J,K)+Lv_CoarsenGlb
    ENDDO
    Lv_OffUniformTot(J,K)=Lv_OffUniform(1,J,K)+&
                          Lv_OffUniform(2,J,K)+Lv_OffUniform(3,J,K)
  ENDDO
ENDDO
min_LvOffUniformTot=MINVAL(Lv_OffUniformTot)
!WRITE(*,21) "MIN LvOffUniformTot in this superblock:",min_LvOffUniformTot
!-----------------------------------------------------------------------!
!#7- Reading index of block matrix of super block
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)
DO K=1,n_BlkSupPre
  READ(10,*)
  READ(10,*)
  READ(10,*) II_Blk(K),JJ_Blk(K),KK_Blk(K)
  WRITE(*,23) "Index of block matrix in super block:",&
              II_Blk(K),JJ_Blk(K),KK_Blk(K)
ENDDO
!-----------------------------------------------------------------------!
!#8- Reading number of overlapping points
READ(10,*)
READ(10,*) n_PointOverlap
!-----------------------------------------------------------------------!
!#9- flag of multi-time-step or not, flag_Mltstp=1: use multi-time-step
READ(10,*)
READ(10,*) flag_Mltstp
!-----------------------------------------------------------------------!
!#10- No. CPUs plan to use, np.
READ(10,*)
READ(10,*) np
!-----------------------------------------------------------------------!
!#11- Ratio of array size between after processing and pre-processing
READ(10,*)
READ(10,*) ratio_PointAftPre
!-----------------------------------------------------------------------!
!#12- lag of more split than initial balance or not, 1 = more split, 0 = no more
READ(10,*)
READ(10,*) flag_MoreSplit
!-----------------------------------------------------------------------!
!#13- No. of extra splits activated by flag_MoreSplit=1.
READ(10,*)
READ(10,*) n_MoreSplit
CLOSE(10)
WRITE(*,'(A50,1X,A2)') "Mesh dimension is:",cas_Dim
WRITE(*,'(A50,1X,A6)') "Mesh type is:",mesh_Type
WRITE(*,*) "Lv_OffUniform(x)	(y)	(z)"
WRITE(*,'(7X,I2,7X,I2,7X,I2)') (((Lv_OffUniform(I,J,K),I=1,3),&
                               J=1,max_BlkGlbPre),K=1,n_BlkSupPre)
!-----------------------------------------------------------------------!
!   Read Mesh Files for BCs in each block: 1st Time     		!
!   Read the 1st time for node info. to allocate			!
!-----------------------------------------------------------------------!
i_BlkGlb=0
DO K=1,n_BlkSupPre
! Mesh type selector
  IF(mesh_Type=="cfx4") THEN
    WRITE(FILENAME,'(A5,I1,A4)') "cfx4-",K,".geo"
    OPEN(20,FILE=FILENAME,STATUS="OLD")
    READ(20,*)
!-----------------------------------------------------------------------!
    tmp_Int=n_BlkPre(K)
    READ(20,*) n_BlkPre(K),n_PathchPre,n_GluePre,n_ElmPre,n_PointPre(K)
    WRITE(message,'(A,I2)') "No. blocks mismatch in superblock:",K
    IF(n_BlkPre(K)/=tmp_Int) CALL PSEXIT(message)
    WRITE(*,22) "Original elements and points:",n_ElmPre,n_PointPre
    max_PointPre=MAXVAL(n_PointPre)
    READ(20,*)
    DO J=1,n_BlkPre(K)
!-----------------------------------------------------------------------!
      ! Add 1 to all node number
      READ(20,*) tmp_Char,II_Pre(J,K),JJ_Pre(J,K),KK_Pre(J,K)
      II_Pre(J,K)=II_Pre(J,K)+1
      JJ_Pre(J,K)=JJ_Pre(J,K)+1
      KK_Pre(J,K)=KK_Pre(J,K)+1
      n_PointLocPre(J,K)=II_Pre(J,K)*JJ_Pre(J,K)*KK_Pre(J,K)
      IF(cas_Dim=="2D") n_PointLocPre(J,K)=n_PointLocPre(J,K)/2
    ENDDO
    CLOSE(20)
  ENDIF
ENDDO
max_PointLocPre=MAXVAL(n_PointLocPre)
ALLOCATE(X_LocPre(max_PointLocPre))
ALLOCATE(Y_LocPre(max_PointLocPre))
ALLOCATE(Z_LocPre(max_PointLocPre))
ALLOCATE(      X_Pre(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(      Y_Pre(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(      Z_Pre(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(     X_toLP(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(     Y_toLP(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(     Z_toLP(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(     X_Leap(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(     Y_Leap(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(     Z_Leap(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(  X_toOvrlp(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(  Y_toOvrlp(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(  Z_toOvrlp(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(  X_Overlap(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(  Y_Overlap(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(  Z_Overlap(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(  X_toSplit(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(  Y_toSplit(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
ALLOCATE(  Z_toSplit(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
! XYZ_Pre = Original, XYZ = after leaping, XYZ_LocPre = Original Local
!-----------------------------------------------------------------------!
!   Read Mesh Files for BCs in each block: 2nd Time     !
!   Read the 2nd time for coordinate.			!
!-----------------------------------------------------------------------!
DO K=1,n_BlkSupPre
! Mesh type selector
  IF(mesh_Type=="cfx4") THEN
    WRITE(FILENAME,'(A5,I1,A4)') "cfx4-",K,".geo"
    OPEN(20,FILE=FILENAME,STATUS="OLD")
    READ(20,*)
!-----------------------------------------------------------------------!
    READ(20,*) n_BlkPre(K),n_PathchPre,n_GluePre,n_ElmPre,n_PointPre(K)
!-----------------------------------------------------------------------!
    READ(20,*)
    DO J=1,n_BlkPre(K)
      READ(20,*)
      ! Already read in last step
    ENDDO
    READ(20,*)
!-----------------------------------------------------------------------!
    DO i_I=1,2*n_PathchPre
      READ(20,*) 
    ENDDO
    READ(20,*)
!-----------------------------------------------------------------------!
    DO i_I=1,n_GluePre
      READ(20,*)
    ENDDO
!-----------------------------------------------------------------------!
    ! Reading mesh, loop on blocks in a super block
    i_Point=0
    i_Point0=0
    i_BlkGlb2=i_BlkGlb
    DO J=1,n_BlkPre(K)
      i_BlkGlb=i_BlkGlb2+J
!      WRITE(*,21) "i_BlkGlb=",i_BlkGlb
!-----------------------------------------------------------------------!
    ! Read X,Y,Z values in each block
      READ(20,*)
      III=II_Pre(J,K)
      JJJ=JJ_Pre(J,K)
      KKK=KK_Pre(J,K)
      DO i_K=1,KKK
        DO i_J=1,JJJ
          DO i_I=1,III
            i_Loc=(i_K-1)*III*JJJ+(i_J-1)*III+i_I
           ! READ(20,*) X_LocPre(i_Loc),Y_LocPre(i_Loc),Z_LocPre(i_Loc)
            i_Point=i_Point0+i_Loc
            READ(20,*) X_Pre(i_Point,K),Y_Pre(i_Point,K),Z_Pre(i_Point,K)
            IF(isidtp_Pre(3,J,K)==3000 .AND. i_J==1) THEN
              Y_Pre(i_Point,K)=0.0044642857 ! This value is half grid length scale
              !Y_LocPre(i_Loc)=Y_LocPre(i_Loc+III)/2 does not work well        
            ENDIF
          !  i_Point=i_Point0+i_Loc
          !  X_Pre(i_Point,K)=X_LocPre(i_Loc)
          !  Y_Pre(i_Point,K)=Y_LocPre(i_Loc)
          !  Z_Pre(i_Point,K)=Z_LocPre(i_Loc)
            !!!!! Half gird processing at axis!!!! Now only for specific case, 2D, axis at bc index 3.
          ENDDO
        ENDDO
      ENDDO
      IF(cas_Dim=="2D") THEN
        i_Point=i_Point-n_PointLocPre(J,K)
        KK_Pre(J,K)=KK_Pre(J,K)/2
      ENDIF
      i_Point0=i_Point
!      WRITE(*,21) "No. points in this block:",i_Point
    ENDDO
  ENDIF
ENDDO
END
