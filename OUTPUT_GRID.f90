!-----------------------------------------------------------------------!
!   Output GRID blocks  		    	 		        !
!   i_BlkGlb: index of block before this output				!
!   i_Point2: index of point before this output				!
!   blk_Suffix: suffix of output file, 6 digit, e.g.:			!
!   _PNTLP,_SPLIT,_OVRLP, started by lower bar "_"
!-----------------------------------------------------------------------!
subroutine output_Grid(&
                       II_toOut,JJ_toOut,KK_toOut,&
                       X_toOut,Y_toOut,Z_toOut,blk_Suffix)
!
USE read_gridcon_pre
USE read_mesh_original
USE system_var
USE point_leaping
USE block_overlaping
!
IMPLICIT NONE
INTEGER I,J,K,IJ_start,IJ_end
INTEGER i_BlkGlb,i_Point2,n_Point
INTEGER II_toOut(max_BlkGlbPre,n_BlkSupPre)
INTEGER JJ_toOut(max_BlkGlbPre,n_BlkSupPre)
INTEGER KK_toOut(max_BlkGlbPre,n_BlkSupPre)
REAL(8) X_toOut(max_PointPre*ratio_PointAftPre,n_BlkSupPre)
REAL(8) Y_toOut(max_PointPre*ratio_PointAftPre,n_BlkSupPre)
REAL(8) Z_toOut(max_PointPre*ratio_PointAftPre,n_BlkSupPre)
REAL(8) X_Out(max_PointPre*ratio_PointAftPre)
REAL(8) Y_Out(max_PointPre*ratio_PointAftPre)
REAL(8) Z_Out(max_PointPre*ratio_PointAftPre)
CHARACTER(LEN=6) blk_Suffix
20 FORMAT(A50)
21 FORMAT(A50,1X,I8)
22 FORMAT(A50,1X,I8,1X,I8)
23 FORMAT(A50,1X,I8,1X,I8,1X,I8)
100 FORMAT(3F16.12)
tmp_Int=max_PointPre*ratio_PointAftPre
i_BlkGlb=0
i_Point2=0
DO K=1,n_BlkSupPre
  DO I=1,tmp_Int
   X_Out(I)=X_toOut(I,K)
   Y_Out(I)=Y_toOut(I,K)
   Z_Out(I)=Z_toOut(I,K)
  ENDDO
  DO J=1,n_BlkPre(K)
    n_Point=II_toOut(J,K)*JJ_toOut(J,K)*KK_toOut(J,K)
    IJ_start=i_Point2+1
    IJ_end=i_Point2+n_Point
    i_BlkGlb=i_BlkGlb+1
    WRITE(FILENAME,'(A7,A6,I3.3)')"GRIDBLO",blk_Suffix,i_BlkGlb
    OPEN(30,FILE=FILENAME,STATUS="UNKNOWN")
    IF(cas_Dim=="2D") THEN
      WRITE(30,*) JJ_toOut(J,K),II_toOut(J,K)
      WRITE(30,100) (X_Out(I),I=IJ_start,IJ_end)
      WRITE(30,100) (Y_Out(I),I=IJ_start,IJ_end)
    ELSE
      WRITE(30,*) KK_toOut(J,K),JJ_toOut(J,K),II_toOut(J,K)
      WRITE(30,100) (X_Out(I),I=IJ_start,IJ_end)
      WRITE(30,100) (Y_Out(I),I=IJ_start,IJ_end)
      WRITE(30,100) (Z_Out(I),I=IJ_start,IJ_end)
    ENDIF
    i_Point2=IJ_end
    CLOSE(30)
  ENDDO
ENDDO
END
