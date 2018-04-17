!-----------------------------------------------------------------------!
!   Perform mesh overlapping		    	 		        !
!-----------------------------------------------------------------------!
subroutine perform_BlkOverlap(i_Run_BlkOvrlp)
!
USE read_gridcon_pre
USE read_mesh_original
USE system_var
USE point_leaping
USE block_overlaping
!
IMPLICIT NONE
! Mesh overlapping
INTEGER I,J,K
INTEGER i_PntOvrlap,i_PntOvrlp,i_PointLoc,i_Swp
INTEGER i_PntLow,i_PntHigh
INTEGER i_I,i_J,i_K,i_l,i_m,I_Leap,J_Leap,K_Leap
INTEGER Diff,Leap, flag_StOvrlp,i_Run_BlkOvrlp
INTEGER gap_WithNxtBlkPreSwpPreSwp,gap_WithNxtPnt,gap_WithNxtPntInNxtBlkPre
INTEGER n_PntThisBlkNow,n_PntBfNxtBlk,Lv_Diff(3),n_PntPreSwp,n_PntAftOvrlpPre
INTEGER II_PreSwp,JJ_PreSwp,KK_PreSwp,n_PointAdd
INTEGER n_PntThisBlkPreSwp
INTEGER II_NxtBlkPreSwp,JJ_NxtBlkPreSwp,KK_NxtBlkPreSwp
INTEGER flag_RunOvrlp,max_IJK_Loc
INTEGER leapX,leapY,leapZ
REAL(8) X_NxtBlkLow,X_NxtBlkHigh,X_YLow,X_YHigh
REAL(8) Y_NxtBlkLow,Y_NxtBlkHigh,Y_YLow,Y_YHigh
REAL(8) Z_NxtBlkLow,Z_NxtBlkHigh,Z_YLow,Z_YHigh
REAL(8) X_PrePntNxtBlkPreSwp,X_NxtPntNxtBlkPreSwp
REAL(8) Y_PrePntNxtBlkPreSwp,Y_NxtPntNxtBlkPreSwp
REAL(8) Z_PrePntNxtBlkPreSwp,Z_NxtPntNxtBlkPreSwp
INTEGER,DIMENSION(:),ALLOCATABLE:: i_PntNxtBlk,nx,ny,nz,nx_l,ny_l,nz_l
INTEGER,DIMENSION(:),ALLOCATABLE:: imx,imy,imz,n_PntAdd
REAL(8),DIMENSION(:),ALLOCATABLE :: X_PreSwp,Y_PreSwp,Z_PreSwp
! n_PointAdd: No. point to be added in one previous step
20 FORMAT(A50)
21 FORMAT(A50,1X,I8)
22 FORMAT(A50,1X,I8,1X,I8)
ALLOCATE(i_PntNxtBlk(n_PointOverlap))
ALLOCATE(n_PntAdd(3))
DO K=1,n_BlkSupPre
  n_PntBfNxtBlk=0
  n_PntAftOvrlpPre=n_PntAftOvrlp(K)
  n_PntAftOvrlp(K)=0
  n_PntPreSwp=0
  i_Swp=0
  tmp_Int=max_PointPre*ratio_PointAftPre   
  ALLOCATE(X_PreSwp(tmp_Int),Y_PreSwp(tmp_Int),Z_PreSwp(tmp_Int))
  ! Assign coordinates after leaping to first overlap coordinate array
  ! X_Overlap(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
!-----------------------------------------------------------------------!
! First overlapping, assign data of leaping to overlap
  IF(i_Run_BlkOvrlp==1) THEN
    tmp_Int=n_PntAftLp(K)
    DO I=1,tmp_Int
      X_Overlap(I,K)=X_Leap(I,K)
      Y_Overlap(I,K)=Y_Leap(I,K)
      Z_Overlap(I,K)=Z_Leap(I,K)
      IF(I==1 .OR. I==tmp_Int) THEN
      WRITE(*,'(A10,2I8,A3,F18.8)') "X_Overlap(",I,K,")=:",X_Overlap(I,K)
      WRITE(*,'(A10,2I8,A3,F18.8)') "Y_Overlap(",I,K,")=:",Y_Overlap(I,K)
      WRITE(*,'(A10,2I8,A3,F18.8)') "Z_Overlap(",I,K,")=:",Z_Overlap(I,K)
      !!! IF(X_Overlap(I,K)/=X_1st .OR. X_Overlap(I,K)/=X_Last) THEN
      !!!   CALL PSEXIT("X_Overlap(1 or n_Pnt,K), first OR last  point error..."
      !!! ELSEIF(Y_Overlap(I,K)/=Y_1st .OR. Y_Overlap(I,K)/=Y_Last) THEN
      !!!   CALL PSEXIT("Y_Overlap(1 or n_Pnt,K), first OR last  point error..."
      !!! ELSEIF(Z_Overlap(I,K)/=Z_1st .OR. Z_Overlap(I,K)/=Z_Last) THEN
      !!!   CALL PSEXIT("Z_Overlap(1 or n_Pnt,K), first OR last  point error..."
      !!! ENDIF
      ENDIF
      X_PreSwp(I)=X_Overlap(I,K)
      Y_PreSwp(I)=Y_Overlap(I,K)
      Z_PreSwp(I)=Z_Overlap(I,K)
    ENDDO
    DO J=1,n_BlkPre(K)
      II_PntOvrlp(J,K)=II_PointLeap(J,K)
      JJ_PntOvrlp(J,K)=JJ_PointLeap(J,K)
      KK_PntOvrlp(J,K)=KK_PointLeap(J,K)
    ENDDO
  ELSE
    ! Not first overlapping, X/Y/Z_Overlap and II/JJ/KK_PntOvrlp(J,K) is known
    DO I=1,n_PntAftOvrlpPre
      X_PreSwp(I)=X_Overlap(I,K)
      Y_PreSwp(I)=Y_Overlap(I,K)
      Z_PreSwp(I)=Z_Overlap(I,K)
    ENDDO
    WRITE(*,*)"X_Overlap(1,1)=:",X_Overlap(1,1)
    WRITE(*,*)"X_Overlap(I,1)=:",X_Overlap(I,1)
  ENDIF
!-----------------------------------------------------------------------!
  DO J=1,n_BlkPre(K)
    max_IJK_Loc=MAX(II_PntOvrlp(J,K),JJ_PntOvrlp(J,K),KK_PntOvrlp(J,K))
    ALLOCATE(nx(max_IJK_Loc),ny(max_IJK_Loc),nz(max_IJK_Loc))
    ALLOCATE(nx_l(max_IJK_Loc),ny_l(max_IJK_Loc),nz_l(max_IJK_Loc))
    ALLOCATE(imx(max_IJK_Loc),imy(max_IJK_Loc),imz(max_IJK_Loc))
    WRITE(*,21)"i_Block_Below=:",J
    flag_RunOvrlp=0
    II_PreSwp=II_PntOvrlp(J,K)
    JJ_PreSwp=JJ_PntOvrlp(J,K)
    KK_PreSwp=KK_PntOvrlp(J,K)
    ! Condition for I/J/K direction stretching to overlap
    IF(i_Run_BlkOvrlp==1) THEN
      gap_WithNxtBlkPreSwpPreSwp=1
    ELSEIF(i_Run_BlkOvrlp==2) THEN
      gap_WithNxtBlkPreSwpPreSwp=II_Blk(K)
    ELSEIF(i_Run_BlkOvrlp==3) THEN
      gap_WithNxtBlkPreSwpPreSwp=II_Blk(K)*JJ_Blk(K)
    ENDIF   
    WRITE(*,21) "gap_WithNxtBlkPreSwpPreSwp=:",gap_WithNxtBlkPreSwpPreSwp
    II_NxtBlkPreSwp=II_PntOvrlp(J+gap_WithNxtBlkPreSwpPreSwp,K)
    JJ_NxtBlkPreSwp=JJ_PntOvrlp(J+gap_WithNxtBlkPreSwpPreSwp,K)
    KK_NxtBlkPreSwp=KK_PntOvrlp(J+gap_WithNxtBlkPreSwpPreSwp,K)
    ! If second boundary surface/edge is fluid domain,then strech layers
    IF(i_Run_BlkOvrlp==1 .AND. isidtp_Pre(2,J,k)==9000) THEN
      ! Condition for I/J/K direction stretching to overlap
      flag_RunOvrlp=1
      II_PntOvrlp(J,K)=II_PntOvrlp(J,K)+n_PointOverlap
      gap_WithNxtPnt=1
      gap_WithNxtPntInNxtBlkPre=1
    ELSEIF(i_Run_BlkOvrlp==2 .AND. isidtp_Pre(4,J,k)==9000) THEN
      flag_RunOvrlp=1
      JJ_PntOvrlp(J,K)=JJ_PntOvrlp(J,K)+n_PointOverlap
      gap_WithNxtPnt=II_PntOvrlp(J,K)
      gap_WithNxtPntInNxtBlkPre=II_NxtBlkPreSwp
    ELSEIF(i_Run_BlkOvrlp==3 .AND. isidtp_Pre(6,J,k)==9000) THEN
      flag_RunOvrlp=1
      KK_PntOvrlp(J,K)=KK_PntOvrlp(J,K)+n_PointOverlap
      gap_WithNxtPnt=II_PntOvrlp(J,K)*JJ_PntOvrlp(J,K)
      gap_WithNxtPntInNxtBlkPre=II_NxtBlkPreSwp*JJ_NxtBlkPreSwp
    ENDIF
    ! n_PntThisBlkPreSwp : total points in this block before streching
    n_PntThisBlkPreSwp=II_PreSwp*JJ_PreSwp*KK_PreSwp
    ! n_PntThisBlkNow : total points in this block after streching
    n_PntThisBlkNow=II_PntOvrlp(J,K)*JJ_PntOvrlp(J,K)*KK_PntOvrlp(J,K)
    IF(n_PntThisBlkNow==n_PntThisBlkPreSwp) THEN
      WRITE(*,*) "--------------------------------------&
                  ----------------------------------------"
      WRITE(*,*) "Block:",J,",Boudnary Blk, No Overlapping Needed!"
      WRITE(*,*) "--------------------------------------&
                  ----------------------------------------"
    ENDIF
    n_PntBfNxtBlk=n_PntBfNxtBlk+n_PntThisBlkPreSwp
!WRITE(*,*) "Lv_OffUniform(x)    (y)     (z)"
!WRITE(*,'(7X,I2,7X,I2,7X,I2)') (((Lv_OffUniform(i_I,i_J,i_K),i_I=1,3),&
!                               i_J=1,max_BlkGlbPre),i_K=1,n_BlkSupPre)
    IF(flag_RunOvrlp==1) THEN
      Lv_Diff(1)=Lv_OffUniform(1,J,K)-Lv_OffUniform(1,J+gap_WithNxtBlkPreSwpPreSwp,K)
      Lv_Diff(2)=Lv_OffUniform(2,J,K)-Lv_OffUniform(2,J+gap_WithNxtBlkPreSwpPreSwp,K)
      Lv_Diff(3)=Lv_OffUniform(3,J,K)-Lv_OffUniform(3,J+gap_WithNxtBlkPreSwpPreSwp,K)
      WRITE(*,*) "lV_diff:",Lv_Diff
      WRITE(*,*)"i/j/k_PreSwp",II_PreSwp,JJ_PreSwp,KK_PreSwp
      IF(II_PreSwp==0) Stop 1
      DO i_K=1,KK_PreSwp
        DO i_J=1,JJ_PreSwp
          DO i_I=1,II_PreSwp
            i_PointLoc=n_PntPreSwp+&
                       (i_K-1)*II_PreSwp*JJ_PreSwp+(i_J-1)*II_PreSwp+i_I
            i_PntOvrlp=n_PntAftOvrlp(K)+&
                       (i_K-1)*II_PntOvrlp(J,K)*JJ_PntOvrlp(J,K)+&
                       (i_J-1)*II_PntOvrlp(J,K)+i_I
            X_Overlap(i_PntOvrlp,K)=X_PreSwp(i_PointLoc)
            Y_Overlap(i_PntOvrlp,K)=Y_PreSwp(i_PointLoc)
            Z_Overlap(i_PntOvrlp,K)=Z_PreSwp(i_PointLoc)
          !  WRITE(*,*)"i_PointLoc,i_PntOvrlp",i_PointLoc,i_PntOvrlp
          !  WRITE(*,*)"X_Overlap(i_PntOvrlp,K)=:",X_Overlap(i_PntOvrlp,K)
            ! Set flag_StOvrlp to 0 before each run.
            flag_StOvrlp=0
            IF(i_I==II_PreSwp) WRITE(*,21) "Last I, i_I=:",i_I
!-----------------------------------------------------------------------!
            IF(i_Run_BlkOvrlp==1 .AND. i_I==II_PreSwp) THEN
              WRITE(*,20) "Reaching I boundary of block, overlapping....."
              WRITE(*,*) "Level difference:I,J,K",Lv_Diff
              flag_StOvrlp=1
              ! Level diff should have the same sign in three directions, all >=0 or all < 0
              IF(Lv_Diff(1)>=0) THEN
                WRITE(*,20) "Coarser than or same with next block, leaping."
                leapX=2**Lv_Diff(1)
                leapY=2**Lv_Diff(2)
                leapZ=2**Lv_Diff(3)
                DO i_l=1,n_PointOverlap
                  ! Calculate point index in next block by leaping, Diff>=0
                  i_PntOvrlp=i_PntOvrlp+gap_WithNxtPnt
                  ! tmp_Int1/2/3 = i/j/k
                  tmp_Int1=i_l*leapX+1
                  tmp_Int2=(i_J-1)*leapY+1
                  tmp_Int3=(i_K-1)*leapZ+1 !! Check Z direction
                  i_PntNxtBlk(i_l)=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                              JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                  X_Overlap(i_PntOvrlp,K)=X_PreSwp(i_PntNxtBlk(i_l))
                  Y_Overlap(i_PntOvrlp,K)=Y_PreSwp(i_PntNxtBlk(i_l))
                  Z_Overlap(i_PntOvrlp,K)=Z_PreSwp(i_PntNxtBlk(i_l))
                  WRITE(*,21)"i_l=:",i_l
                  WRITE(*,21)"i_PntOvrlp=:",i_PntOvrlp
                  WRITE(*,21)"i_PntNxtBlk=:",i_PntNxtBlk(i_l)
                  WRITE(*,*)"X=:",X_Overlap(i_PntOvrlp,K)
                  WRITE(*,*)"Y=:",Y_Overlap(i_PntOvrlp,K)
                  WRITE(*,*)"Z=:",Z_Overlap(i_PntOvrlp,K)
                ENDDO
              ELSEIF(Lv_Diff(1)<0) THEN
                WRITE(*,20) "Finer than next block, overlap by adding."
                n_PntAdd(1)=2**((-1)*Lv_Diff(1))-1
                n_PntAdd(2)=2**((-1)*Lv_Diff(2))-1
                n_PntAdd(3)=2**((-1)*Lv_Diff(3))-1
                DO i_l=1,n_PointOverlap
                  ! Calculate point index in next block by finding the integer index
                  ! nx,ny,nz = Integer index in adding point
                  ! nx/y/z_l = Lower integer of nx/y/z if nx/y/z is not integer
                  ! imx/y/z  = remainder
                  i_PntOvrlp=i_PntOvrlp+gap_WithNxtPnt
                  ! Swip I, x
                  nx(i_l)=i_l/(n_PntAdd(1)+1)+1
                  IF((nx(i_l)-1)*(n_PntAdd(1)+1)==i_l) THEN
                    nx_l(i_l)=0
                    imx(i_l)=0
                  ELSE
                    nx_l(i_l)=nx(i_l)
                    nx(i_l)=0
                    imx(i_l)=i_l-(nx_l(i_l)-1)*(n_PntAdd(1)+1)
                  ENDIF
                  ! Swip I, y
                  ny(i_J)=(i_J-1)/(n_PntAdd(2)+1)+1
                  IF((ny(i_J)-1)*(n_PntAdd(2)+1)+1==i_J) THEN
                    ny_l(i_J)=0
                    imy(i_J)=0
                  ELSE
                    ny_l(i_J)=ny(i_J)
                    ny(i_J)=0
                    imy(i_J)=i_J-1-(ny_l(i_J)-1)*(n_PntAdd(2)+1)
                  ENDIF
                  ! Swip I, z
                  nz(i_K)=(i_K-1)/(n_PntAdd(3)+1)+1
                  IF((nz(i_K)-1)*(n_PntAdd(3)+1)+1==i_K) THEN
                    nz_l(i_K)=0
                    imz(i_K)=0
                  ELSE
                    nz_l(i_K)=nz(i_K)
                    nz(i_K)=0
                    imz(i_K)=i_K-1-(nz_l(i_K)-1)*(n_PntAdd(3)+1)
                  ENDIF
                  IF(imx(i_l)==0 .AND. imy(i_J)==0) THEN ! 2D case, Swipe I
                    tmp_Int1=nx(i_l)
                    tmp_Int2=ny(i_J)
                    tmp_Int3=i_K
                    i_PntNxtBlk(i_l)=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    X_Overlap(i_PntOvrlp,K)=X_PreSwp(i_PntNxtBlk(i_l))
                    Y_Overlap(i_PntOvrlp,K)=Y_PreSwp(i_PntNxtBlk(i_l))
                    Z_Overlap(i_PntOvrlp,K)=Z_PreSwp(i_PntNxtBlk(i_l))
                    WRITE(*,*) "X_Overlap(",i_PntOvrlp,")",X_Overlap(i_PntOvrlp,K)
                  ELSEIF(imx(i_l)==0 .AND. imy(i_J)/=0) THEN
                    tmp_Int1=nx(i_l)
                    tmp_Int2=ny_l(i_J)
                    tmp_Int3=i_K
                    i_PntLow=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    tmp_Int2=ny_l(i_J)+1
                    i_PntHigh=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    X_NxtBlkLow=X_PreSwp(i_PntLow)
                    Y_NxtBlkLow=Y_PreSwp(i_PntLow)
                    Z_NxtBlkLow=Z_PreSwp(i_PntLow)
                    X_NxtBlkHigh=X_PreSwp(i_PntHigh)
                    Y_NxtBlkHigh=Y_PreSwp(i_PntHigh)
                    Z_NxtBlkHigh=Z_PreSwp(i_PntHigh)
                    X_Overlap(i_PntOvrlp,K)=X_NxtBlkLow+imy(i_J)*&
                                (X_NxtBlkHigh-X_NxtBlkLow)/(n_PntAdd(2)+1)
                    Y_Overlap(i_PntOvrlp,K)=Y_NxtBlkLow+imy(i_J)*&
                                (Y_NxtBlkHigh-Y_NxtBlkLow)/(n_PntAdd(2)+1)
                    Z_Overlap(i_PntOvrlp,K)=Z_NxtBlkLow+imy(i_J)*&
                                (Z_NxtBlkHigh-Z_NxtBlkLow)/(n_PntAdd(2)+1)
                  ELSEIF(imx(i_l)/=0 .AND. imy(i_J)==0) THEN
                    tmp_Int1=nx_l(i_l)
                    tmp_Int2=ny(i_J)
                    tmp_Int3=i_K
                    i_PntLow=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    tmp_Int1=nx_l(i_l)+1
                    i_PntHigh=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    X_NxtBlkLow=X_PreSwp(i_PntLow)
                    Y_NxtBlkLow=Y_PreSwp(i_PntLow)
                    Z_NxtBlkLow=Z_PreSwp(i_PntLow)
                    X_NxtBlkHigh=X_PreSwp(i_PntHigh)
                    Y_NxtBlkHigh=Y_PreSwp(i_PntHigh)
                    Z_NxtBlkHigh=Z_PreSwp(i_PntHigh)
                    X_Overlap(i_PntOvrlp,K)=X_NxtBlkLow+imx(i_l)*&
                                (X_NxtBlkHigh-X_NxtBlkLow)/(n_PntAdd(1)+1)
                    Y_Overlap(i_PntOvrlp,K)=Y_NxtBlkLow+imx(i_l)*&
                                (Y_NxtBlkHigh-Y_NxtBlkLow)/(n_PntAdd(1)+1)
                    Z_Overlap(i_PntOvrlp,K)=Z_NxtBlkLow+imx(i_l)*&
                                (Z_NxtBlkHigh-Z_NxtBlkLow)/(n_PntAdd(1)+1)
                  ELSEIF(imx(i_l)/=0 .AND. imy(i_J)/=0) THEN
                    ! X-Y- AND X+Y-
                    tmp_Int1=nx_l(i_l)
                    tmp_Int2=ny_l(i_J)
                    tmp_Int3=i_K
                    i_PntLow=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    tmp_Int1=nx_l(i_l)+1
                    i_PntHigh=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    X_NxtBlkLow=X_PreSwp(i_PntLow)
                    Y_NxtBlkLow=Y_PreSwp(i_PntLow)
                    Z_NxtBlkLow=Z_PreSwp(i_PntLow)
                    X_NxtBlkHigh=X_PreSwp(i_PntHigh)
                    Y_NxtBlkHigh=Y_PreSwp(i_PntHigh)
                    Z_NxtBlkHigh=Z_PreSwp(i_PntHigh)
                    X_YLow=X_NxtBlkLow+imx(i_l)*&
                                (X_NxtBlkHigh-X_NxtBlkLow)/(n_PntAdd(1)+1)
                    Y_YLow=Y_NxtBlkLow+imx(i_l)*&
                                (Y_NxtBlkHigh-Y_NxtBlkLow)/(n_PntAdd(1)+1)
                    Z_YLow=Z_NxtBlkLow+imx(i_l)*&
                                (Z_NxtBlkHigh-Z_NxtBlkLow)/(n_PntAdd(1)+1)
                    ! X-Y+ AND X+Y+
                    tmp_Int1=nx_l(i_l)
                    tmp_Int2=ny_l(i_J)+1
                    tmp_Int3=i_K
                    i_PntLow=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    tmp_Int1=nx_l(i_l)+1
                    i_PntHigh=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    X_NxtBlkLow=X_PreSwp(i_PntLow)
                    Y_NxtBlkLow=Y_PreSwp(i_PntLow)
                    Z_NxtBlkLow=Z_PreSwp(i_PntLow)
                    X_NxtBlkHigh=X_PreSwp(i_PntHigh)
                    Y_NxtBlkHigh=Y_PreSwp(i_PntHigh)
                    Z_NxtBlkHigh=Z_PreSwp(i_PntHigh)
                    X_YHigh=X_NxtBlkLow+imx(i_l)*&
                                (X_NxtBlkHigh-X_NxtBlkLow)/(n_PntAdd(1)+1)
                    Y_YHigh=Y_NxtBlkLow+imx(i_l)*&
                                (Y_NxtBlkHigh-Y_NxtBlkLow)/(n_PntAdd(1)+1)
                    Z_YHigh=Z_NxtBlkLow+imx(i_l)*&
                                (Z_NxtBlkHigh-Z_NxtBlkLow)/(n_PntAdd(1)+1)
                    X_Overlap(i_PntOvrlp,K)=X_YLow+imy(i_J)*&
                                (X_YHigh-X_YLow)/(n_PntAdd(2)+1)
                    Y_Overlap(i_PntOvrlp,K)=Y_YLow+imy(i_J)*&
                                (Y_YHigh-Y_YLow)/(n_PntAdd(2)+1)
                    Z_Overlap(i_PntOvrlp,K)=Z_YLow+imy(i_J)*&
                                (Z_YHigh-Z_YLow)/(n_PntAdd(2)+1)
                  ENDIF
                ENDDO
              ENDIF
!-----------------------------------------------------------------------!
            ELSEIF(i_Run_BlkOvrlp==2 .AND. i_J==JJ_PreSwp) THEN
              WRITE(*,20) "Reaching J boundary of block, overlapping....."
              WRITE(*,21) "Level difference:",Lv_Diff(2)
              flag_StOvrlp=1
              ! Level diff should have the same sign in three directions, all >=0 or all < 0
              IF(Lv_Diff(2)>=0) THEN
                WRITE(*,20) "Coarser than or same with next block, leaping."
                leapX=2**Lv_Diff(1)
                leapY=2**Lv_Diff(2)
                leapZ=2**Lv_Diff(3)
                DO i_l=1,n_PointOverlap
                  ! Calculate point index in next block by leaping, Diff>=0
                  i_PntOvrlp=i_PntOvrlp+gap_WithNxtPnt
                  ! tmp_Int1/2/3 = i/j/k
                  tmp_Int1=(i_I-1)*leapX+1
                  tmp_Int2=i_l*leapY+1
                  tmp_Int3=(i_K-1)*leapZ+1 !! Check Z direction
                  i_PntNxtBlk(i_l)=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                              JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                  X_Overlap(i_PntOvrlp,K)=X_PreSwp(i_PntNxtBlk(i_l))
                  Y_Overlap(i_PntOvrlp,K)=Y_PreSwp(i_PntNxtBlk(i_l))
                  Z_Overlap(i_PntOvrlp,K)=Z_PreSwp(i_PntNxtBlk(i_l))
                ENDDO
              ELSEIF(Lv_Diff(2)<0) THEN
                WRITE(*,20) "Finer than next block, overlap by adding."
                n_PntAdd(1)=2**((-1)*Lv_Diff(1))-1
                n_PntAdd(2)=2**((-1)*Lv_Diff(2))-1
                n_PntAdd(3)=2**((-1)*Lv_Diff(3))-1
                DO i_l=1,n_PointOverlap
                  ! Calculate point index in next block by finding the integer index
                  i_PntOvrlp=i_PntOvrlp+gap_WithNxtPnt
                  ! Swip J, x
                  nx(i_I)=(i_I-1)/(n_PntAdd(1)+1)+1
                  IF((nx(i_I)-1)*(n_PntAdd(1)+1)+1==i_I) THEN
                    nx_l(i_I)=0
                    imx(i_I)=0
                  ELSE
                    nx_l(i_I)=nx(i_I)
                    nx(i_I)=0
                    imx(i_I)=i_I-1-(nx_l(i_I)-1)*(n_PntAdd(1)+1)
                  ENDIF
                  ! Swip J, y
                  ny(i_l)=i_l/(n_PntAdd(2)+1)+1
                  IF((ny(i_l)-1)*(n_PntAdd(2)+1)==i_l) THEN
                    ny_l(i_l)=0
                    imy(i_l)=0
                  ELSE
                    ny_l(i_l)=ny(i_l)
                    ny(i_l)=0
                    imy(i_l)=i_l-(ny_l(i_l)-1)*(n_PntAdd(2)+1)
                  ENDIF
                  ! Swip J, z
                  nz(i_K)=(i_K-1)/(n_PntAdd(3)+1)+1
                  IF((nz(i_K)-1)*(n_PntAdd(3)+1)+1==i_K) THEN
                    nz_l(i_K)=0
                    imz(i_K)=0
                  ELSE
                    nz_l(i_K)=nz(i_K)
                    nz(i_K)=0
                    imz(i_K)=i_K-1-(nz_l(i_K)-1)*(n_PntAdd(3)+1)
                  ENDIF
                  IF(imx(i_I)==0 .AND. imy(i_l)==0) THEN ! 2D case, Swipe J
                    tmp_Int1=nx(i_I)
                    tmp_Int2=ny(i_l)
                    tmp_Int3=i_K
                    i_PntNxtBlk(i_l)=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    X_Overlap(i_PntOvrlp,K)=X_PreSwp(i_PntNxtBlk(i_l))
                    Y_Overlap(i_PntOvrlp,K)=Y_PreSwp(i_PntNxtBlk(i_l))
                    Z_Overlap(i_PntOvrlp,K)=Z_PreSwp(i_PntNxtBlk(i_l))
                  ELSEIF(imx(i_I)==0 .AND. imy(i_l)/=0) THEN
                    tmp_Int1=nx(i_I)
                    tmp_Int2=ny_l(i_l)
                    tmp_Int3=i_K
                    i_PntLow=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    tmp_Int2=ny_l(i_l)+1
                    i_PntHigh=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    X_NxtBlkLow=X_PreSwp(i_PntLow)
                    Y_NxtBlkLow=Y_PreSwp(i_PntLow)
                    Z_NxtBlkLow=Z_PreSwp(i_PntLow)
                    X_NxtBlkHigh=X_PreSwp(i_PntHigh)
                    Y_NxtBlkHigh=Y_PreSwp(i_PntHigh)
                    Z_NxtBlkHigh=Z_PreSwp(i_PntHigh)
                    X_Overlap(i_PntOvrlp,K)=X_NxtBlkLow+imy(i_l)*&
                                (X_NxtBlkHigh-X_NxtBlkLow)/(n_PntAdd(2)+1)
                    Y_Overlap(i_PntOvrlp,K)=Y_NxtBlkLow+imy(i_l)*&
                                (Y_NxtBlkHigh-Y_NxtBlkLow)/(n_PntAdd(2)+1)
                    Z_Overlap(i_PntOvrlp,K)=Z_NxtBlkLow+imy(i_l)*&
                                (Z_NxtBlkHigh-Z_NxtBlkLow)/(n_PntAdd(2)+1)
                  ELSEIF(imx(i_I)/=0 .AND. imy(i_l)==0) THEN
                    tmp_Int1=nx_l(i_I)
                    tmp_Int2=ny(i_l)
                    tmp_Int3=i_K
                    i_PntLow=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    tmp_Int1=nx_l(i_I)+1
                    i_PntHigh=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    X_NxtBlkLow=X_PreSwp(i_PntLow)
                    Y_NxtBlkLow=Y_PreSwp(i_PntLow)
                    Z_NxtBlkLow=Z_PreSwp(i_PntLow)
                    X_NxtBlkHigh=X_PreSwp(i_PntHigh)
                    Y_NxtBlkHigh=Y_PreSwp(i_PntHigh)
                    Z_NxtBlkHigh=Z_PreSwp(i_PntHigh)
                    X_Overlap(i_PntOvrlp,K)=X_NxtBlkLow+imx(i_I)*&
                                (X_NxtBlkHigh-X_NxtBlkLow)/(n_PntAdd(1)+1)
                    Y_Overlap(i_PntOvrlp,K)=Y_NxtBlkLow+imx(i_I)*&
                                (Y_NxtBlkHigh-Y_NxtBlkLow)/(n_PntAdd(1)+1)
                    Z_Overlap(i_PntOvrlp,K)=Z_NxtBlkLow+imx(i_I)*&
                                (Z_NxtBlkHigh-Z_NxtBlkLow)/(n_PntAdd(1)+1)
                  ELSEIF(imx(i_I)/=0 .AND. imy(i_l)/=0) THEN
                    ! X-Y- AND X+Y-
                    tmp_Int1=nx_l(i_I)
                    tmp_Int2=ny_l(i_l)
                    tmp_Int3=i_K
                    i_PntLow=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    tmp_Int1=nx_l(i_I)+1
                    i_PntHigh=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    X_NxtBlkLow=X_PreSwp(i_PntLow)
                    Y_NxtBlkLow=Y_PreSwp(i_PntLow)
                    Z_NxtBlkLow=Z_PreSwp(i_PntLow)
                    X_NxtBlkHigh=X_PreSwp(i_PntHigh)
                    Y_NxtBlkHigh=Y_PreSwp(i_PntHigh)
                    Z_NxtBlkHigh=Z_PreSwp(i_PntHigh)
                    X_YLow=X_NxtBlkLow+imx(i_I)*&
                                (X_NxtBlkHigh-X_NxtBlkLow)/(n_PntAdd(1)+1)
                    Y_YLow=Y_NxtBlkLow+imx(i_I)*&
                                (Y_NxtBlkHigh-Y_NxtBlkLow)/(n_PntAdd(1)+1)
                    Z_YLow=Z_NxtBlkLow+imx(i_I)*&
                                (Z_NxtBlkHigh-Z_NxtBlkLow)/(n_PntAdd(1)+1)
                    ! X-Y+ AND X+Y+
                    tmp_Int1=nx_l(i_I)
                    tmp_Int2=ny_l(i_J)+1
                    tmp_Int3=i_K
                    i_PntLow=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    tmp_Int1=nx_l(i_I)+1
                    i_PntHigh=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                                JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                    X_NxtBlkLow=X_PreSwp(i_PntLow)
                    Y_NxtBlkLow=Y_PreSwp(i_PntLow)
                    Z_NxtBlkLow=Z_PreSwp(i_PntLow)
                    X_NxtBlkHigh=X_PreSwp(i_PntHigh)
                    Y_NxtBlkHigh=Y_PreSwp(i_PntHigh)
                    Z_NxtBlkHigh=Z_PreSwp(i_PntHigh)
                    X_YHigh=X_NxtBlkLow+imx(i_I)*&
                                (X_NxtBlkHigh-X_NxtBlkLow)/(n_PntAdd(1)+1)
                    Y_YHigh=Y_NxtBlkLow+imx(i_I)*&
                                (Y_NxtBlkHigh-Y_NxtBlkLow)/(n_PntAdd(1)+1)
                    Z_YHigh=Z_NxtBlkLow+imx(i_I)*&
                                (Z_NxtBlkHigh-Z_NxtBlkLow)/(n_PntAdd(1)+1)
                    X_Overlap(i_PntOvrlp,K)=X_YLow+imy(i_l)*&
                                (X_YHigh-X_YLow)/(n_PntAdd(2)+1)
                    Y_Overlap(i_PntOvrlp,K)=Y_YLow+imy(i_l)*&
                                (Y_YHigh-Y_YLow)/(n_PntAdd(2)+1)
                    Z_Overlap(i_PntOvrlp,K)=Z_YLow+imy(i_l)*&
                                (Z_YHigh-Z_YLow)/(n_PntAdd(2)+1)
                  ENDIF
                ENDDO
              ENDIF
!-----------------------------------------------------------------------!
            ELSEIF(i_Run_BlkOvrlp==3 .AND. i_K==KK_PreSwp) THEN
              WRITE(*,20) "Reaching K boundary of block, overlapping....."
              WRITE(*,21) "Level difference:",Lv_Diff(3)
              flag_StOvrlp=1
              ! Level diff should have the same sign in three directions, all >=0 or all < 0
              IF(Lv_Diff(3)>=0) THEN
                WRITE(*,20) "Coarser than or same with next block, leaping."
                leapX=2**Lv_Diff(1)
                leapY=2**Lv_Diff(2)
                leapZ=2**Lv_Diff(3)
                DO i_l=1,n_PointOverlap
                  ! Calculate point index in next block by leaping, Diff>=0
                  i_PntOvrlp=i_PntOvrlp+gap_WithNxtPnt
                  ! tmp_Int1/2/3 = i/j/k
                  tmp_Int1=(i_I-1)*leapX+1
                  tmp_Int2=(i_J-1)*leapY+1
                  tmp_Int3=i_l*leapZ+1 !! Check Z direction
                  i_PntNxtBlk(i_l)=n_PntBfNxtBlk+(tmp_Int3-1)*II_NxtBlkPreSwp*&
                              JJ_NxtBlkPreSwp+(tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                  X_Overlap(i_PntOvrlp,K)=X_PreSwp(i_PntNxtBlk(i_l))
                  Y_Overlap(i_PntOvrlp,K)=Y_PreSwp(i_PntNxtBlk(i_l))
                  Z_Overlap(i_PntOvrlp,K)=Z_PreSwp(i_PntNxtBlk(i_l))
                ENDDO
              ELSEIF(Lv_Diff(3)<0) THEN
                ! This module, to be developed
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ELSE
!-----------------------------------------------------------------------!
      ! If no overlapping is to be performed in some blocks.
      ! e.g. boundary blocks and corner blocks.
      WRITE(*,21)"No. points before i_Block:",J
      WRITE(*,21)"In PreSwp:",n_PntPreSwp
      WRITE(*,21)"In AftSwp:",n_PntAftOvrlp(K)
      DO i_K=1,KK_PreSwp
        DO i_J=1,JJ_PreSwp
          DO i_I=1,II_PreSwp
            tmp_Int1=(i_K-1)*II_PreSwp*JJ_PreSwp+(i_J-1)*II_PreSwp+i_I
            tmp_Int2=(i_K-1)*II_PntOvrlp(J,K)*JJ_PntOvrlp(J,K)+&
                     (i_J-1)*II_PntOvrlp(J,K)+i_I
            IF(tmp_Int1/=tmp_Int2) THEN
              CALL PSEXIT("Error in point overlapping when flag_RunOvrlp=0")
            ENDIF
            i_PointLoc=n_PntPreSwp+tmp_Int1
            i_PntOvrlp=n_PntAftOvrlp(K)+tmp_Int2
            X_Overlap(i_PntOvrlp,K)=X_PreSwp(i_PointLoc)
            Y_Overlap(i_PntOvrlp,K)=Y_PreSwp(i_PointLoc)
            Z_Overlap(i_PntOvrlp,K)=Z_PreSwp(i_PointLoc)
          ENDDO
        ENDDO
      ENDDO
    ENDIF    
    n_PntPreSwp=n_PntPreSwp+n_PntThisBlkPreSwp
    n_PntAftOvrlp(K)=n_PntAftOvrlp(K)+n_PntThisBlkNow
    WRITE(*,21)"No. points after i_Block:",J
    WRITE(*,21)"In PreSwp:",n_PntPreSwp
    WRITE(*,21)"In AftSwp:",n_PntAftOvrlp(K)
    DEALLOCATE(nx,ny,nz)
    DEALLOCATE(nx_l,ny_l,nz_l)
    DEALLOCATE(imx,imy,imz)
  ENDDO
ENDDO
END
