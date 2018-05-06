!-----------------------------------------------------------------------!
!   Perform block splitting		    	 		        !
!   Split block into 2**n_Split sub-blocks by one2two method		!
!-----------------------------------------------------------------------!
subroutine perform_BlkSplit
!
USE read_gridcon_pre
USE read_mesh_original
USE system_var
USE point_leaping
USE block_overlaping
USE block_splitting
!
IMPLICIT NONE
! Block splitting
INTEGER I,J,K,i_I,i_J,i_K,i_Point,i_PointPre
INTEGER min_OffUniformInSup,n_PntEqvGlbPreSwp,avg_PntPerCpuPreSwp
INTEGER n_BlkInBlkPreSwp,n_BlkGlbAftSplit
INTEGER i_BlkGlb,i_BlkGlbPre,i_Point2
INTEGER II,JJ,KK,IJ_start,IJ_end,isidtp_Aft(6),Lv_OffUnfrmAft(3)
INTEGER i_BlkSplit,i_Split
INTEGER i_BlkSplit1,i_BlkSplit2
INTEGER i_BlkSplit11,i_BlkSplit12,i_BlkSplit21,i_BlkSplit22
INTEGER i_BlkSplit111,i_BlkSplit112,i_BlkSplit121,i_BlkSplit122,&
        i_BlkSplit211,i_BlkSplit212,i_BlkSplit221,i_BlkSplit222
INTEGER i_BlkSplit1111,i_BlkSplit1112,i_BlkSplit1121,i_BlkSplit1122,&
        i_BlkSplit1211,i_BlkSplit1212,i_BlkSplit1221,i_BlkSplit1222,&
        i_BlkSplit2111,i_BlkSplit2112,i_BlkSplit2121,i_BlkSplit2122,&
        i_BlkSplit2211,i_BlkSplit2212,i_BlkSplit2221,i_BlkSplit2222
INTEGER i_BlkSplit11111,i_BlkSplit11112,i_BlkSplit11121,i_BlkSplit11122,&
        i_BlkSplit11211,i_BlkSplit11212,i_BlkSplit11221,i_BlkSplit11222,&
        i_BlkSplit12111,i_BlkSplit12112,i_BlkSplit12121,i_BlkSplit12122,&
        i_BlkSplit12211,i_BlkSplit12212,i_BlkSplit12221,i_BlkSplit12222,&
        i_BlkSplit21111,i_BlkSplit21112,i_BlkSplit21121,i_BlkSplit21122,&
        i_BlkSplit21211,i_BlkSplit21212,i_BlkSplit21221,i_BlkSplit21222,&
        i_BlkSplit22111,i_BlkSplit22112,i_BlkSplit22121,i_BlkSplit22122,&
        i_BlkSplit22211,i_BlkSplit22212,i_BlkSplit22221,i_BlkSplit22222
REAL(8) dx,dy,dz
INTEGER,DIMENSION(:),ALLOCATABLE :: n_PntEqvPreSwp
INTEGER,DIMENSION(:,:),ALLOCATABLE :: np_LocInt,np_Loc2nInt
REAL(8),DIMENSION(:),ALLOCATABLE :: X_Out,Y_Out,Z_Out
REAL(8),DIMENSION(:,:),ALLOCATABLE :: np_LocReal

20 FORMAT(A50)
21 FORMAT(A50,1X,I8)
22 FORMAT(A50,1X,I8,1X,I8)
23 FORMAT(A50,1X,I8,1X,I8,1X,I8)
100   FORMAT(2F16.12)
ALLOCATE(                    min_dxyz(n_BlkSupPre))
ALLOCATE(   min_dxyzLoc(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(                      min_dt(n_BlkSupPre))
ALLOCATE(            dt(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(        n_Step(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(   n_PntEqvLoc(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(              n_PntEqvPreSwp(n_BlkSupPre))
ALLOCATE(     np_LocInt(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(    np_LocReal(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(       n_Split(max_BlkGlbPre,n_BlkSupPre))
ALLOCATE(   np_Loc2nInt(max_BlkGlbPre,n_BlkSupPre))
!-----------------------------------------------------------------------!
! Calculate the mininmum edge in each superblock
min_dxyz=100.0
min_dxyzLoc=100.0
DO K=1,n_BlkSupPre
  i_PointPre=0
  DO J=1,n_BlkPre(K)
    !FIND THE MIN Dxyz in each superblock where has the same min_dxyz
    DO i_K=1,KK_toSplit(J,K)
      DO i_J=1,JJ_toSplit(J,K)
        DO i_I=1,II_toSplit(J,K)
          i_Point=i_PointPre+(i_K-1)*II_toSplit(J,K)*JJ_toSplit(J,K)+&
                  (i_J-1)*II_toSplit(J,K)+i_I       
          IF(i_I<II_toSplit(J,K)) THEN
            dx=DABS(X_toSplit(i_Point+1,K)-X_toSplit(i_Point,K))
          ELSE
            dx=100.0
          ENDIF
          IF(i_J<JJ_toSplit(J,K)) THEN          
            dy=DABS(Y_toSplit(i_Point+II_toSplit(J,K),K)-Y_toSplit(i_Point,K))
          ELSE
            dy=100.0
          ENDIF
          IF(i_K<KK_toSplit(J,K)) THEN          
            dz=DABS(Z_toSplit(i_Point+II_toSplit(J,K)*JJ_toSplit(J,K),K)-&
                  Z_toSplit(i_Point,K)) 
          ELSE
            dz=100.0
          ENDIF
          tmp_Real=MIN(dx,dy,dz)
          IF(cas_Dim=="2D") tmp_Real=MIN(dx,dy)
          IF(tmp_Real<min_dxyzLoc(J,K)) THEN
            min_dxyzLoc(J,K)=tmp_Real
          ENDIF
          IF(tmp_Real<min_dxyz(K)) THEN
            min_dxyz(K)=tmp_Real
          ENDIF
          IF(min_dxyz(K)<1e-10) THEN
            WRITE(*,21)"min_dxyz = 0 in i_Blk=:",J
            WRITE(*,21)"I=:",i_I
            WRITE(*,21)"J=:",i_J
            WRITE(*,21)"K=:",i_K
            WRITE(*,*) "i_Point:",i_Point
            WRITE(*,*) "min_dxyz:",min_dxyz(K)
            CALL PSEXIT("min_dxyz got value 0!")
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    i_PointPre=i_Point
  ENDDO
ENDDO
!-----------------------------------------------------------------------!
!Calculate dt and number of steps to march on for a unit time.
CFL=1.D2
n_Step=0
min_OffUniformInSup=100
! Check if all elements become 0
DO K=1,n_BlkSupPre
  min_dt(K)=CFL*min_dxyz(K)
  DO J=1,n_BlkPre(K)
        IF(flag_Mltstp==1) THEN
          dt(J,K)=CFL*min_dxyzLoc(J,K)
          WRITE(*,*) "Multi-time-step: dt = ", dt(J,K)
        ELSE
          dt(J,K)=min_dt(K)
          IF(J==1) WRITE(*,*) "Single-time-step: dt = ", dt(J,K)
        ENDIF
        n_Step(J,K)=1/dt(J,K)
  ENDDO
ENDDO
max_N_Step=MAXVAL(n_step)
!WRITE(*,21)"max_N_Step:",max_N_Step
!-----------------------------------------------------------------------!
!Calculate average equivalent points per CPU before swiping
n_PntEqvLoc=0
n_PntEqvPreSwp=0
n_PntEqvGlbPreSwp=0
DO K=1,n_BlkSupPre
  DO J=1,n_BlkPre(K)
    n_PntInBlk(J,K)=II_toSplit(J,K)*JJ_toSplit(J,K)*KK_toSplit(J,K)
    n_PntEqvLoc(J,K)=n_PntInBlk(J,K)*REAL(n_Step(J,K))/REAL(max_N_Step)
    WRITE(*,200) "i_Block:",J,"n_PntEqvLoc(J,K):",n_PntEqvLoc(J,K)
200 FORMAT(A13,I8.8,A29,I8)
    n_PntEqvPreSwp(K)=n_PntEqvPreSwp(K)+n_PntEqvLoc(J,K)
  ENDDO
  n_PntEqvGlbPreSwp=n_PntEqvGlbPreSwp+n_PntEqvPreSwp(K)
ENDDO
avg_PntPerCpuPreSwp=n_PntEqvGlbPreSwp/np
WRITE(*,21)"avg_PntPerCpuPreSwP(Equivalent):",avg_PntPerCpuPreSwp
!-----------------------------------------------------------------------!
! Calculate total blocks after block splitting and do array allocation
n_BlkGlbAftSplit=0
DO K=1,n_BlkSupPre
  DO J=1,n_BlkPre(K)
    n_Split(J,K)=0
    np_LocReal(J,K)=REAL(n_PntEqvLoc(J,K))/REAL(avg_PntPerCpuPreSwp)
    np_LocInt(J,K)=np_LocReal(J,K)
    n_BlkGlbAftSplit=n_BlkGlbAftSplit+1
!    WRITE(*,*) "np_LocReal,np_LocInt:",np_LocReal(J,K),np_LocInt(J,K)
    IF(np_LocReal(J,K)>1) THEN
      ! This block is larger than 1.5 times of average, perform splitting
      IF(np_LocReal(J,K)-np_LocInt(J,K)>=0.5 .AND. np_LocInt(J,K)<2) THEN
        np_LocInt(J,K)=np_LocInt(J,K)+1
      ELSEIF(np_LocReal(J,K)-np_LocInt(J,K)>=0.8 .AND. np_LocInt(J,K)<7) THEN
        np_LocInt(J,K)=np_LocInt(J,K)+1
      ENDIF
      DO WHILE(2**n_Split(J,K)<np_LocInt(J,K))
        n_Split(J,K)=n_Split(J,K)+1
      ENDDO
      IF(n_Split(J,K)>=4) THEN
      IF(2**n_Split(J,K)-np_LocInt(J,K) >          &
         1.5*(2**(n_Split(J,K)-1)-np_LocInt(J,K))) &
         n_Split(J,K)=n_Split(J,K)-1
      ENDIF
      np_Loc2nInt(J,K)=2**n_Split(J,K)
      n_BlkGlbAftSplit=n_BlkGlbAftSplit+np_Loc2nInt(J,K)-1
    ENDIF
    IF(flag_MoreSplit==1) THEN
      n_Split(J,K)=n_Split(J,K)+n_MoreSplit
      IF(n_Split(J,K)>6) &
        CALL PSEXIT("No. Split larger than 6, decrease n_MoreSplit.")
    ENDIF
  ENDDO
ENDDO
WRITE(*,20) "n_Split:"
WRITE(*,'(I2,1X,I2,1X,I2,1X,I2)'),n_Split
!-----------------------------------------------------------------------!
! Perform block splitting
i_BlkGlb=0
i_BlkGlbPre=0
i_Point2=0
DO K=1,n_BlkSupPre
  WRITE(*,*) "--------------------------------------&
            ----------------------------------------"
  WRITE(*,21)"Process block splitting, i_SupBlk=:",K
  DO J=1,n_BlkPre(K)
    n_toSplitLoc(J,K)=II_toSplit(J,K)*JJ_toSplit(J,K)*KK_toSplit(J,K)
    i_BlkGlbPre=i_BlkGlbPre+1
    WRITE(*,*) "--------------------------------------&
            ----------------------------------------"
    WRITE(*,21)"Process block splitting, i_Block=:",J
!    WRITE(*,21)"i_BlockGlb=:",i_BlkGlbPre
!    WRITE(*,21)"n_Split=:",n_Split(J,K)
!-----------------------------------------------------------------------!
!   Export blocks for which no splitting is needed
    IJ_start=i_Point2+1
    IJ_end=i_Point2+n_toSplitLoc(J,K)
    i_Point2=IJ_end
    IF(n_Split(J,K)==0) THEN
      II=II_toSplit(J,K)
      JJ=JJ_toSplit(J,K)
      KK=KK_toSplit(J,K)
      ALLOCATE(X_Out(n_toSplitLoc(J,K)))
      ALLOCATE(Y_Out(n_toSplitLoc(J,K)))
      ALLOCATE(Z_Out(n_toSplitLoc(J,K)))
      X_Out(1:n_toSplitLoc(J,K))=X_toSplit(IJ_start:IJ_end,K)
      Y_Out(1:n_toSplitLoc(J,K))=Y_toSplit(IJ_start:IJ_end,K)
      Z_Out(1:n_toSplitLoc(J,K))=Z_toSplit(IJ_start:IJ_end,K)
      i_BlkGlb=i_BlkGlb+1
      isidtp_Aft(1:6)=isidtp_Pre(1:6,J,K)
      Lv_OffUnfrmAft(1:3)=Lv_OffUniform(1:3,J,K)
      CALL function_Output_Grid(II,JJ,KK,i_BlkGlb,&
            X_Out,Y_Out,Z_Out,n_toSplitLoc(J,K),  &
            isidtp_Aft,Lv_OffUnfrmAft,dt(J,K),n_PntEqvLoc(J,K))
      DEALLOCATE(X_Out,Y_Out,Z_Out)
    ENDIF
!-----------------------------------------------------------------------!
!   Start to split block needed to be split
    DO i_Split=1,n_Split(J,K)
      !Loop on splitting cuts.
      SELECT CASE(i_Split)
      CASE(1)
        ! Split 1
        i_BlkSplit=i_BlkGlbPre
        CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
      CASE(2)
      ! Split 2
        tmp_Int=i_BlkSplit
        i_BlkSplit1=tmp_Int*10+1
        i_BlkSplit=i_BlkSplit1
!        WRITE(*,*)"Lv_2-1:"
        CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit2=tmp_Int*10+2
        i_BlkSplit=i_BlkSplit2
!        WRITE(*,*)"Lv_2-2:"
        CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
      CASE(3) ! One into 8 blocks
        ! 11,12,21,22
        i_BlkSplit11=i_BlkSplit1*10+1
        i_BlkSplit=i_BlkSplit11
!        WRITE(*,*)"i_BlkSplit:",i_BlkSplit 
!        WRITE(*,*)"Lv_3-1:"
        CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit12=i_BlkSplit1*10+2
        i_BlkSplit=i_BlkSplit12
!        WRITE(*,*)"i_BlkSplit:",i_BlkSplit 
!        WRITE(*,*)"Lv_3-2:"
        CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit21=i_BlkSplit2*10+1
        i_BlkSplit=i_BlkSplit21
!        WRITE(*,*)"i_BlkSplit:",i_BlkSplit 
!        WRITE(*,*)"Lv_3-3:"
        CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit22=i_BlkSplit2*10+2
        i_BlkSplit=i_BlkSplit22
!        WRITE(*,*)"i_BlkSplit:",i_BlkSplit 
!        WRITE(*,*)"Lv_3-4:"
        CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
      CASE(4) ! One into 16 blocks
        ! 111,112,121,122,211,212,221,222
        i_BlkSplit111=i_BlkSplit11*10+1
        i_BlkSplit=i_BlkSplit111
401     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit112=i_BlkSplit11*10+2
        i_BlkSplit=i_BlkSplit112
402     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit121=i_BlkSplit12*10+1
        i_BlkSplit=i_BlkSplit121
403     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit122=i_BlkSplit12*10+2
        i_BlkSplit=i_BlkSplit122
404     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit211=i_BlkSplit21*10+1
        i_BlkSplit=i_BlkSplit211
405     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit212=i_BlkSplit21*10+2
        i_BlkSplit=i_BlkSplit212
406     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit221=i_BlkSplit22*10+1
        i_BlkSplit=i_BlkSplit221
407      CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit222=i_BlkSplit22*10+2
        i_BlkSplit=i_BlkSplit222
408     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
      CASE(5) !One into 32 blocks
        !1111,1112
        i_BlkSplit1111=i_BlkSplit111*10+1
        i_BlkSplit=i_BlkSplit1111
501     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit1112=i_BlkSplit111*10+2
        i_BlkSplit=i_BlkSplit1112
502     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !1121,1122
        i_BlkSplit1121=i_BlkSplit112*10+1
        i_BlkSplit=i_BlkSplit1121
503     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit1122=i_BlkSplit112*10+2
        i_BlkSplit=i_BlkSplit1122
504     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !1211,1212 
        i_BlkSplit1211=i_BlkSplit121*10+1
        i_BlkSplit=i_BlkSplit1211
505     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit1212=i_BlkSplit121*10+2
        i_BlkSplit=i_BlkSplit1212
506     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !1221,1222
        i_BlkSplit1221=i_BlkSplit122*10+1
        i_BlkSplit=i_BlkSplit1221
507     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit1222=i_BlkSplit122*10+2
        i_BlkSplit=i_BlkSplit1222
508     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !2111,2112
        i_BlkSplit2111=i_BlkSplit211*10+1
        i_BlkSplit=i_BlkSplit2111
509     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit2112=i_BlkSplit211*10+2
        i_BlkSplit=i_BlkSplit2112
510     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !2121,2122
        i_BlkSplit2121=i_BlkSplit212*10+1
        i_BlkSplit=i_BlkSplit2121
511     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit2122=i_BlkSplit212*10+2
        i_BlkSplit=i_BlkSplit2122
512     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        ! 2211,2212
        i_BlkSplit2211=i_BlkSplit221*10+1
        i_BlkSplit=i_BlkSplit2211
513     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit2212=i_BlkSplit221*10+2
        i_BlkSplit=i_BlkSplit2212
514     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        ! 2221,2222
        i_BlkSplit2221=i_BlkSplit222*10+1
        i_BlkSplit=i_BlkSplit2221
515     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit2222=i_BlkSplit222*10+2
        i_BlkSplit=i_BlkSplit2222
516     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
      CASE(6) ! Split one block into 64 blocks
        !11111,11112
        i_BlkSplit11111=i_BlkSplit1111*10+1
        i_BlkSplit=i_BlkSplit11111
601     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit11112=i_BlkSplit1111*10+2
        i_BlkSplit=i_BlkSplit11112
602     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !11121,11122
        i_BlkSplit11121=i_BlkSplit1112*10+1
        i_BlkSplit=i_BlkSplit11121
603     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit11122=i_BlkSplit1112*10+2
        i_BlkSplit=i_BlkSplit11122
604     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !11211,11212
        i_BlkSplit11211=i_BlkSplit1121*10+1
        i_BlkSplit=i_BlkSplit11211
605     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit11212=i_BlkSplit1121*10+2
        i_BlkSplit=i_BlkSplit11212
606     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !11221,11222
        i_BlkSplit11221=i_BlkSplit1122*10+1
        i_BlkSplit=i_BlkSplit11221
607     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit11222=i_BlkSplit1122*10+2
        i_BlkSplit=i_BlkSplit11222
608     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !12111,12112
        i_BlkSplit12111=i_BlkSplit1211*10+1
        i_BlkSplit=i_BlkSplit12111
609     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit12112=i_BlkSplit1211*10+2
        i_BlkSplit=i_BlkSplit12112
610     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !12121,12122
        i_BlkSplit12121=i_BlkSplit1212*10+1
        i_BlkSplit=i_BlkSplit12121
611     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit12122=i_BlkSplit1212*10+2
        i_BlkSplit=i_BlkSplit12122
612     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !12211,12212
        i_BlkSplit12211=i_BlkSplit1221*10+1
        i_BlkSplit=i_BlkSplit12211
613     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit12212=i_BlkSplit1221*10+2
        i_BlkSplit=i_BlkSplit12212
614     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !12221,12222
        i_BlkSplit12221=i_BlkSplit1222*10+1
        i_BlkSplit=i_BlkSplit12221
615     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit12222=i_BlkSplit1222*10+2
        i_BlkSplit=i_BlkSplit12222
616     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !21111,21112
        i_BlkSplit21111=i_BlkSplit2111*10+1
        i_BlkSplit=i_BlkSplit21111
617     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit21112=i_BlkSplit2111*10+2
        i_BlkSplit=i_BlkSplit21112
618     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !21121,21122
        i_BlkSplit21121=i_BlkSplit2112*10+1
        i_BlkSplit=i_BlkSplit21121
619     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit21122=i_BlkSplit2112*10+2
        i_BlkSplit=i_BlkSplit21122
620     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !21211,21212
        i_BlkSplit21211=i_BlkSplit2121*10+1
        i_BlkSplit=i_BlkSplit21211
621     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit21212=i_BlkSplit2121*10+2
        i_BlkSplit=i_BlkSplit21212
622     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !21221,21222
        i_BlkSplit21221=i_BlkSplit2122*10+1
        i_BlkSplit=i_BlkSplit21221
623     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit21222=i_BlkSplit2122*10+2
        i_BlkSplit=i_BlkSplit21222
624     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !22111,22112
        i_BlkSplit22111=i_BlkSplit2211*10+1
        i_BlkSplit=i_BlkSplit22111
625     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit22112=i_BlkSplit2211*10+2
        i_BlkSplit=i_BlkSplit22112
626     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !22121,22122
        i_BlkSplit22121=i_BlkSplit2212*10+1
        i_BlkSplit=i_BlkSplit22121
627     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit22122=i_BlkSplit2212*10+2
        i_BlkSplit=i_BlkSplit22122
628     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !22211,22212
        i_BlkSplit22211=i_BlkSplit2221*10+1
        i_BlkSplit=i_BlkSplit22211
629     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit22212=i_BlkSplit2221*10+2
        i_BlkSplit=i_BlkSplit22212
630     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        !22221,22222
        i_BlkSplit22221=i_BlkSplit2222*10+1
        i_BlkSplit=i_BlkSplit22221
631     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
        i_BlkSplit22222=i_BlkSplit2222*10+2
        i_BlkSplit=i_BlkSplit22222
632     CALL function_BlkSplit(i_BlkSplit,i_Split,i_BlkGlb,J,K,IJ_start,IJ_end)
      CASE DEFAULT
        WRITE(*,*) "No. split is not enough, please increase it!"
        CALL PSEXIT("No. split is not enough, please increase it!")
      ENDSELECT
    ENDDO
  ENDDO
ENDDO
n_BlkFinal=i_BlkGlb
WRITE(*,21) "Final No. block, n_BlkFinal:",n_BlkFinal
END
