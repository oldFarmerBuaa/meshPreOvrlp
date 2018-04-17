!-----------------------------------------------------------------------!
!    Perform point leaping						!
!-----------------------------------------------------------------------!
subroutine perform_PntLeap
!
USE read_gridcon_pre
USE read_mesh_original
USE system_var
USE point_leaping
USE block_overlaping
!
INTEGER I,J,K
INTEGER i_Point,i_PointPre,i_PointLP,i_PointLoc
INTEGER I_Leap,J_Leap,K_Leap
INTEGER i_Loc,i_I_aft,i_J_aft,i_K_aft,i_Loop,j_Loop,k_Loop,i_l,n_PointAdd
INTEGER step_PntLP(3)
REAL(8),DIMENSION(:),ALLOCATABLE   :: X_Loc,Y_Loc,Z_Loc

20 FORMAT(A50)
21 FORMAT(A50,1X,I8)
22 FORMAT(A50,1X,I8,1X,I8)
23 FORMAT(A50,1X,I8,1X,I8,1X,I8)
100   FORMAT(2F16.12)
DO K=1,n_BlkSupPre
  WRITE(*,*) "--------------------------------------&
            ----------------------------------------"
  WRITE(*,21) "Process point leaping, i_SupBlk=:",K
  i_PointPre=0
  i_PointLP=0
  i_Point=0
  DO J=1,n_BlkPre(K)
    WRITE(*,*) "--------------------------------------&
            ----------------------------------------"
    WRITE(*,21) "Process point leaping, i_Block=:",J
    ! Calculate leaping step and node number after point leaping
    ! If Lv_OffUniform(1,J,K)>=0, mesh to be coarsened, decrease points
    ! March on with several steps
    ! Otherwise, mesh to be refined, march on with step = 1
    DO I=1,3
      IF(Lv_OffUniform(I,J,K)>=0) THEN
        step_PntLP(I)=2**Lv_OffUniform(I,J,K)
      ELSEIF(Lv_OffUniform(I,J,K)<0) THEN
        step_PntLP(I)=1
      ENDIF
    ENDDO
!    WRITE(*,*) "Lv_OffUniform:",(Lv_OffUniform(I,J,K),I=1,3)
    IF(Lv_OffUniform(1,J,K)>=0) THEN
      II_PntLP(J,K)=(II_toLP(J,K)-1)/step_PntLP(1)+1
    ELSE
      II_PntLP(J,K)=(II_toLP(J,K)-1)*2**((-1)*Lv_OffUniform(1,J,K))+1
    ENDIF
    IF(Lv_OffUniform(2,J,K)>0) THEN
      JJ_PntLP(J,K)=(JJ_toLP(J,K)-1)/step_PntLP(2)+1
    ELSE
      JJ_PntLP(J,K)=(JJ_toLP(J,K)-1)*2**((-1)*Lv_OffUniform(2,J,K))+1
    ENDIF
    IF(cas_Dim=="2D") THEN
      KK_PntLP(J,K)=1
    ELSEIF(Lv_OffUniform(3,J,K)>0) THEN
      KK_PntLP(J,K)=(KK_toLP(J,K)-1)/step_PntLP(3)+1
    ELSE
      KK_PntLP(J,K)=(KK_toLP(J,K)-1)*2**((-1)*Lv_OffUniform(3,J,K))+1
    ENDIF
!    WRITE(*,23) "I,J,K in this block before leaping:",&
!                II_toLP(J,K),JJ_toLP(J,K),KK_toLP(J,K)
    n_PntInBlk(J,K)=II_toLP(J,K)*JJ_toLP(J,K)*KK_toLP(J,K)
!    WRITE(*,23) "I,J,K in this block after leaping:",&
!                II_PntLP(J,K),JJ_PntLP(J,K),KK_PntLP(J,K)
    n_PointLP(J,K)=II_PntLP(J,K)*JJ_PntLP(J,K)*KK_PntLP(J,K)
    IF(n_PointLP(J,K)>5*n_PointPre(K)) THEN
      WRITE(*,*) "Point overflow warning in X,Y,Z after leaping, &
                 increase ratio_PointAftPre"
      STOP 1
    ENDIF
    I_Leap=step_PntLP(1)
    J_Leap=step_PntLP(2)
    K_Leap=step_PntLP(3)
    tmp_Int=n_PointLP(J,K)
    WRITE(*,21) "Total Points in this block before leaping:",&
                n_PntInBlk(J,K)
    WRITE(*,21) "Total Points in this block after leaping:",tmp_Int
!    WRITE(*,23) "Leaping steps:",I_Leap,J_Leap,K_Leap
    ! Total point in each block, 2 layers counted for 2D cases
    ALLOCATE(X_Loc(tmp_Int))
    ALLOCATE(Y_Loc(tmp_Int))
    ALLOCATE(Z_Loc(tmp_Int))
    tmp_Int=0
    k_Loop=1
    DO i_K=1,KK_toLP(J,K),K_Leap
      j_Loop=1
      DO i_J=1,JJ_toLP(J,K),J_Leap
        i_Loop=1
        DO i_I=1,II_toLP(J,K),I_Leap
          i_PointPre=i_Point+(i_K-1)*II_toLP(J,K)*JJ_toLP(J,K)+&
                     (i_J-1)*II_toLP(J,K)+i_I
          IF(Lv_OffUniform(1,J,K)<0) THEN
          i_I_aft=1+(i_Loop-1)*2**((-1)*Lv_OffUniform(1,J,K))
          ELSE
          i_I_aft=i_Loop
          ENDIF
          IF(Lv_OffUniform(2,J,K)<0) THEN
          i_J_aft=1+(j_Loop-1)*2**((-1)*Lv_OffUniform(2,J,K))
          ELSE
          i_J_aft=j_Loop
          ENDIF
          IF(Lv_OffUniform(3,J,K)<0) THEN
          i_K_aft=1+(k_Loop-1)*2**((-1)*Lv_OffUniform(3,J,K))
          ELSE
          i_k_aft=k_Loop
          ENDIF
          i_PointLoc=(i_K_aft-1)*II_PntLP(J,K)*JJ_PntLP(J,K)+&
                     (i_J_aft-1)*II_PntLP(J,K)+i_I_aft
!          WRITE(*,22)"i_PointPre,i_PointLoc:",i_PointPre,i_PointLoc
          X_Loc(i_PointLoc)=X_toLP(i_PointPre,K)
          Y_Loc(i_PointLoc)=Y_toLP(i_PointPre,K)
          Z_Loc(i_PointLoc)=Z_toLP(i_PointPre,K)
          ! Calculate to add inner points
          ! #1---- Add to I grid line, J,K grid line is to be developed.
          IF(Lv_OffUniform(1,J,K)<0) THEN
            IF(i_J_aft==tmp_Int) THEN
            ! If on the same grid line, J = constant
            n_PointAdd=2**((-1)*Lv_OffUniform(1,J,K))-1
            ! n_PointAdd = number of inner points
              DO i_l=1,n_PointAdd
                X_Loc(i_PointLoc-i_l)=X_Loc(i_PointLoc)-&
                                  i_l*(X_Loc(i_PointLoc)-tmp_Real1)/&
                                  (n_PointAdd+1)
                Y_Loc(i_PointLoc-i_l)=Y_Loc(i_PointLoc)-&
                                  i_l*(Y_Loc(i_PointLoc)-tmp_Real2)/&
                                  (n_PointAdd+1)
                Z_Loc(i_PointLoc-i_l)=Z_Loc(i_PointLoc)-&
                                  i_l*(Z_Loc(i_PointLoc)-tmp_Real3)/&
                                  (n_PointAdd+1)
              ENDDO
            ENDIF
          ENDIF
          i_Loop=i_Loop+1
          tmp_Int=i_J_aft
          tmp_Real1=X_Loc(i_PointLoc)
          tmp_Real2=Y_Loc(i_PointLoc)
          tmp_Real3=Z_Loc(i_PointLoc)
          ! For index I use, only
        ENDDO
        j_Loop=j_Loop+1
      ENDDO
      k_Loop=k_Loop+1
    ENDDO
    i_Point2=i_PointLP
    DO i_K=1,KK_PntLP(J,K)
      DO i_J=1,JJ_PntLP(J,K)
        DO i_I=1,II_PntLP(J,K)
          i_PointLoc=(i_K-1)*II_PntLP(J,K)*JJ_PntLP(J,K)+&
                     (i_J-1)*II_PntLP(J,K)+i_I
          i_PointLP=i_Point2+i_PointLoc
          X_Leap(i_PointLP,K)=X_Loc(i_PointLoc)
          Y_Leap(i_PointLP,K)=Y_Loc(i_PointLoc)
          Z_Leap(i_PointLP,K)=Z_Loc(i_PointLoc)
        ENDDO
      ENDDO
    ENDDO
    i_Point=i_Point+n_PntInBlk(J,K)
    DEALLOCATE(X_Loc,Y_Loc,Z_Loc)
  ENDDO
  n_PntAftLp(K)=i_PointLP
  IF(n_PntAftLp(K)>max_PointPre*ratio_PointAftPre) THEN
    WRITE(*,20) "No. point overflew after Leaping,increase ratio_PointAftPre!!!"
    STOP 1
  ENDIF
  WRITE(*,21) "Total points in this super block n_PntAftLp=:",n_PntAftLp(K)
ENDDO
END
