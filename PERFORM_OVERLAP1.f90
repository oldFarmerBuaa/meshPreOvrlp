!-----------------------------------------------------------------------!
!   Perform mesh overlapping		    	 		        !
!-----------------------------------------------------------------------!
subroutine perform_BlkOverlap1(i_RBO)
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
INTEGER i_PntOvrlap,i_PntOvrlp,i_PointLoc
INTEGER i_PntLow,i_PntHigh
INTEGER i_I,i_J,i_K,i_l,i_m,I_Leap,J_Leap,K_Leap
INTEGER Diff,Leap, flag_StOvrlp,i_RBO
INTEGER gap_WithNxtBlkPreSwpPreSwp,gap_WithNxtPnt,gap_WithNxtPntInNxtBlkPre
INTEGER n_PntThisBlkNow,Lv_Diff(3),n_PntPreSwp,n_PntAftOvrlpPre
INTEGER II_PreSwp,JJ_PreSwp,KK_PreSwp,n_PointAdd
INTEGER n_PntThisBlkPreSwp
INTEGER II_NxtBlkPreSwp,JJ_NxtBlkPreSwp,KK_NxtBlkPreSwp
INTEGER flag_RunOvrlp,max_IJK_Loc
INTEGER leapX,leapY,leapZ
INTEGER i_BlkNxt
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
! i_RBO, index of running of Blk overlap
20 FORMAT(A50)
21 FORMAT(A50,1X,I8)
22 FORMAT(A50,1X,I8,1X,I8)
ALLOCATE(i_PntNxtBlk(n_PointOverlap))
ALLOCATE(n_PntAdd(3))
!!!! For overlapping first cases Lv_Diff=0, Uniform mesh needed.
Lv_Diff(1:3)=0
DO K=1,n_BlkSupPre
  WRITE(*,*) "--------------------------------------&
            ----------------------------------------"
  WRITE(*,21)"Process block overlapping, i_SupBlock=:",K
  n_PntAftOvrlpPre=n_PntAftOvrlp(K)
  n_PntAftOvrlp(K)=0
  n_PntPreSwp=0
  tmp_Int=max_PointPre*ratio_PointAftPre   
  ALLOCATE(X_PreSwp(tmp_Int),Y_PreSwp(tmp_Int),Z_PreSwp(tmp_Int))
  ! Assign coordinates after leaping to first overlap coordinate array
  ! X_Overlap(max_PointPre*ratio_PointAftPre,n_BlkSupPre))
!-----------------------------------------------------------------------!
! First overlapping, assign data of leaping to overlap
  IF(i_RBO==1) THEN
    tmp_Int=n_PntToOvrlp(K)
    DO I=1,tmp_Int
      X_Overlap(I,K)=X_toOvrlp(I,K)
      Y_Overlap(I,K)=Y_toOvrlp(I,K)
      Z_Overlap(I,K)=Z_toOvrlp(I,K)
      IF(I==1 .OR. I==tmp_Int) THEN
!      WRITE(*,'(A10,2I8,A3,F18.8)') "X_Overlap(",I,K,")=:",X_Overlap(I,K)
!      WRITE(*,'(A10,2I8,A3,F18.8)') "Y_Overlap(",I,K,")=:",Y_Overlap(I,K)
!      WRITE(*,'(A10,2I8,A3,F18.8)') "Z_Overlap(",I,K,")=:",Z_Overlap(I,K)
      ENDIF
      X_PreSwp(I)=X_Overlap(I,K)
      Y_PreSwp(I)=Y_Overlap(I,K)
      Z_PreSwp(I)=Z_Overlap(I,K)
    ENDDO
    DO J=1,n_BlkPre(K)
      II_PntOvrlp(J,K)=II_toOvrlp(J,K)
      JJ_PntOvrlp(J,K)=JJ_toOvrlp(J,K)
      KK_PntOvrlp(J,K)=KK_toOvrlp(J,K)
    ENDDO
  ELSE
    ! Not first overlapping, X/Y/Z_Overlap and II/JJ/KK_PntOvrlp(J,K) is known
    DO I=1,n_PntAftOvrlpPre
      X_PreSwp(I)=X_Overlap(I,K)
      Y_PreSwp(I)=Y_Overlap(I,K)
      Z_PreSwp(I)=Z_Overlap(I,K)
    ENDDO
  ENDIF
  DO J=1,n_BlkPre(K)
    n_PntInBlk(J,K)=II_PntOvrlp(J,K)*JJ_PntOvrlp(J,K)*KK_PntOvrlp(J,K)
  ENDDO
  n_PntBfBlk=0
  DO J=1,n_BlkPre(K)-1
    n_PntBfBlk(J+1,K)=n_PntBfBlk(J,K)+n_PntInBLK(J,K)
  ENDDO
!-----------------------------------------------------------------------!
  DO J=1,n_BlkPre(K)
    max_IJK_Loc=MAX(II_PntOvrlp(J,K),JJ_PntOvrlp(J,K),KK_PntOvrlp(J,K))
    ALLOCATE(nx(max_IJK_Loc),ny(max_IJK_Loc),nz(max_IJK_Loc))
    ALLOCATE(nx_l(max_IJK_Loc),ny_l(max_IJK_Loc),nz_l(max_IJK_Loc))
    ALLOCATE(imx(max_IJK_Loc),imy(max_IJK_Loc),imz(max_IJK_Loc))
    WRITE(*,*) "--------------------------------------&
            ----------------------------------------"
    WRITE(*,21)"Process block overlapping, i_Block=:",J
    flag_RunOvrlp=0
    II_PreSwp=II_PntOvrlp(J,K)
    JJ_PreSwp=JJ_PntOvrlp(J,K)
    KK_PreSwp=KK_PntOvrlp(J,K)
    ! n_PntThisBlkPreSwp : total points in this block before streching
    n_PntThisBlkPreSwp=II_PreSwp*JJ_PreSwp*KK_PreSwp
    ! Condition for I/J/K direction stretching to overlap
    IF(i_RBO==1) THEN
      gap_WithNxtBlkPreSwpPreSwp=1
    ELSEIF(i_RBO==2) THEN
      gap_WithNxtBlkPreSwpPreSwp=II_Blk(K)
    ELSEIF(i_RBO==3) THEN
      gap_WithNxtBlkPreSwpPreSwp=II_Blk(K)*JJ_Blk(K)
    ENDIF   
    i_BlkNxt=J+gap_WithNxtBlkPreSwpPreSwp
    II_NxtBlkPreSwp=II_PntOvrlp(i_BlkNxt,K)
    JJ_NxtBlkPreSwp=JJ_PntOvrlp(i_BlkNxt,K)
    KK_NxtBlkPreSwp=KK_PntOvrlp(i_BlkNxt,K)
    ! If second boundary surface/edge is fluid domain,then strech layers
    IF(i_RBO==1 .AND. isidtp_Pre(2,J,k)==9000) THEN
      ! Condition for I/J/K direction stretching to overlap
      flag_RunOvrlp=1
      II_PntOvrlp(J,K)=II_PntOvrlp(J,K)+n_PointOverlap
      gap_WithNxtPnt=1
      gap_WithNxtPntInNxtBlkPre=1
    ELSEIF(i_RBO==2 .AND. isidtp_Pre(4,J,k)==9000) THEN
      flag_RunOvrlp=1
      JJ_PntOvrlp(J,K)=JJ_PntOvrlp(J,K)+n_PointOverlap
      gap_WithNxtPnt=II_PntOvrlp(J,K)
      gap_WithNxtPntInNxtBlkPre=II_NxtBlkPreSwp
    ELSEIF(i_RBO==3 .AND. isidtp_Pre(6,J,k)==9000) THEN
      flag_RunOvrlp=1
      KK_PntOvrlp(J,K)=KK_PntOvrlp(J,K)+n_PointOverlap
      gap_WithNxtPnt=II_PntOvrlp(J,K)*JJ_PntOvrlp(J,K)
      gap_WithNxtPntInNxtBlkPre=II_NxtBlkPreSwp*JJ_NxtBlkPreSwp
    ENDIF
    ! n_PntThisBlkNow : total points in this block after streching
    n_PntThisBlkNow=II_PntOvrlp(J,K)*JJ_PntOvrlp(J,K)*KK_PntOvrlp(J,K)
    IF(n_PntThisBlkNow==n_PntThisBlkPreSwp) THEN
!      WRITE(*,*) "--------------------------------------&
!                  ----------------------------------------"
!      WRITE(*,*) "Block:",J,",Boudnary Blk, No Overlapping Needed!"
!      WRITE(*,*) "--------------------------------------&
!                  ----------------------------------------"
    ENDIF
    IF(flag_RunOvrlp==1) THEN
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
            IF(i_RBO==1 .AND. i_I==II_PreSwp) THEN
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
                i_PntNxtBlk(i_l)=n_PntBfBlk(i_BlkNxt,K)+&
                                 (tmp_Int3-1)*II_NxtBlkPreSwp*JJ_NxtBlkPreSwp+&
                                 (tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                X_Overlap(i_PntOvrlp,K)=X_PreSwp(i_PntNxtBlk(i_l))
                Y_Overlap(i_PntOvrlp,K)=Y_PreSwp(i_PntNxtBlk(i_l))
                Z_Overlap(i_PntOvrlp,K)=Z_PreSwp(i_PntNxtBlk(i_l))
              ENDDO
!-----------------------------------------------------------------------!
            ELSEIF(i_RBO==2 .AND. i_J==JJ_PreSwp) THEN
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
                i_PntNxtBlk(i_l)=n_PntBfBlk(i_BlkNxt,K)+&
                                 (tmp_Int3-1)*II_NxtBlkPreSwp*JJ_NxtBlkPreSwp+&
                                 (tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                X_Overlap(i_PntOvrlp,K)=X_PreSwp(i_PntNxtBlk(i_l))
                Y_Overlap(i_PntOvrlp,K)=Y_PreSwp(i_PntNxtBlk(i_l))
                Z_Overlap(i_PntOvrlp,K)=Z_PreSwp(i_PntNxtBlk(i_l))
              ENDDO
!-----------------------------------------------------------------------!
            ELSEIF(i_RBO==3 .AND. i_K==KK_PreSwp) THEN
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
                i_PntNxtBlk(i_l)=n_PntBfBlk(i_BlkNxt,K)+&
                                 (tmp_Int3-1)*II_NxtBlkPreSwp*JJ_NxtBlkPreSwp+&
                                 (tmp_Int2-1)*II_NxtBlkPreSwp+tmp_Int1
                X_Overlap(i_PntOvrlp,K)=X_PreSwp(i_PntNxtBlk(i_l))
                Y_Overlap(i_PntOvrlp,K)=Y_PreSwp(i_PntNxtBlk(i_l))
                Z_Overlap(i_PntOvrlp,K)=Z_PreSwp(i_PntNxtBlk(i_l))
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ELSE
!-----------------------------------------------------------------------!
      ! If no overlapping is to be performed in some blocks.
      ! e.g. boundary blocks and corner blocks.
 !     WRITE(*,21)"No. points before i_Block:",J
 !     WRITE(*,21)"In PreSwp:",n_PntPreSwp
 !     WRITE(*,21)"In AftSwp:",n_PntAftOvrlp(K)
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
!    WRITE(*,21)"No. points after i_Block:",J
!    WRITE(*,21)"In PreSwp:",n_PntPreSwp
!    WRITE(*,21)"In AftSwp:",n_PntAftOvrlp(K)
    DEALLOCATE(nx,ny,nz)
    DEALLOCATE(nx_l,ny_l,nz_l)
    DEALLOCATE(imx,imy,imz)
  ENDDO
  IF(n_PntAftOvrlp(K)>max_PointPre*ratio_PointAftPre) THEN
    WRITE(*,20) "No. point overflew after overlaping,increase ratio_PointAftPre!!!"
    STOP 1
  ENDIF
ENDDO
END
