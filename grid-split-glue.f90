program grid_split_glue
!-----------------------------------------------------------------------!
!	Read, split and glue mesh from ICEM				!
!	1) Same block mesh to multiblock.			 	!
!	2) Output mesh as cfx-4 in file cfx4-N.geo			!
!	3) Run this program to split and glue mesh.			!
!									!
!	Loc =  In a block						!
!	    =  In a super block						!
!	Glb =  In the whole mesh					!
!	Xiaolong Tang 01/28/2018					!
!									!
!	Run sequence I: Overlap -> PointLeap -> Split			!
!		    II: PointLeap -> Overlap -> Split (N/A right now)	!
!	Pre-condition I: Original I/J/K index 				!
!	                 shall meet I/J/K+n_PntOvrlp = 2^N+1		!
!		     II: Original I/J/K index shall meet I/J/K=2^N+1	!
!		         N = No. Grid Levels				!
!-----------------------------------------------------------------------!
!
USE read_gridcon_pre
USE read_mesh_original
USE system_var
USE point_leaping
USE block_overlaping
USE block_splitting
!
IMPLICIT NONE

20 FORMAT(A80)
WRITE(*,20) "=============================== Processing Started ========&
             ====================="
WRITE(*,20) "Input overlap pattern, 1 = One-Direction, 2 = Two-Direction&
             ....................."
READ(*,*) f_Ovrlp1Or2Dir
CALL SYSTEM("rm GRIDBLO*")
CALL SYSTEM("rm output_GRID/GRIDBLO*")
!CALL SYSTEM("if [ ! -d output_GRID ]; mkdir output_GRID")
!-----------------------------------------------------------------------!
!    Read mesh								!
!-----------------------------------------------------------------------!
WRITE(*,20) "===========================================================&
             ====================="
WRITE(*,20) "Start to read mesh files...................................&
             ....................."
CALL read_GrdAndAllocate
WRITE(*,20) "Reading mesh files.........................................&
             .................Done"
CALL output_Grid(II_Pre,JJ_Pre,KK_Pre,&
                 X_Pre,Y_Pre,Z_Pre,"_INPUT")
!-----------------------------------------------------------------------!
!   Perform mesh overlapping --- 1	    	 		        !
!-----------------------------------------------------------------------!
WRITE(*,20) "===========================================================&
             ====================="
WRITE(*,20) "Start to overlap...........................................&
             ....................."
II_toOvrlp=II_Pre
JJ_toOvrlp=JJ_Pre
KK_toOvrlp=KK_Pre
X_toOvrlp=X_Pre
Y_toOvrlp=Y_Pre
Z_toOvrlp=Z_Pre
n_PntToOvrlp=n_PointPre
IF (f_Ovrlp1Or2Dir==1) THEN
CALL perform_BlkOverlap1(1)
CALL output_Grid(II_PntOvrlp,JJ_PntOvrlp,KK_PntOvrlp,&
                 X_Overlap,Y_Overlap,Z_Overlap,"_OVRLP")
WRITE(*,20) "Overlap in I direction.....................................&
             .................Done"
CALL perform_BlkOverlap1(2)
CALL output_Grid(II_PntOvrlp,JJ_PntOvrlp,KK_PntOvrlp,&
                 X_Overlap,Y_Overlap,Z_Overlap,"_OVRLP")
WRITE(*,20) "Overlap in J direction.....................................&
             .................Done"
IF(cas_Dim=="3D") THEN
  CALL perform_BlkOverlap1(3)
  CALL output_Grid(II_PntOvrlp,JJ_PntOvrlp,KK_PntOvrlp,&
                 X_Overlap,Y_Overlap,Z_Overlap,"_OVRLP")
WRITE(*,20) "Overlap in K direction.....................................&
             .................Done"
ENDIF
ELSEIF (f_Ovrlp1Or2Dir==2) THEN
CALL perform_BlkOverlap2(1)
CALL output_Grid(II_PntOvrlp,JJ_PntOvrlp,KK_PntOvrlp,&
                 X_Overlap,Y_Overlap,Z_Overlap,"_OVRLP")
WRITE(*,20) "Overlap in I direction.....................................&
             .................Done"
CALL perform_BlkOverlap2(2)
CALL output_Grid(II_PntOvrlp,JJ_PntOvrlp,KK_PntOvrlp,&
                 X_Overlap,Y_Overlap,Z_Overlap,"_OVRLP")
WRITE(*,20) "Overlap in J direction.....................................&
             .................Done"
IF(cas_Dim=="3D") THEN
  CALL perform_BlkOverlap2(3)
  CALL output_Grid(II_PntOvrlp,JJ_PntOvrlp,KK_PntOvrlp,&
                 X_Overlap,Y_Overlap,Z_Overlap,"_OVRLP")
WRITE(*,20) "Overlap in K direction.....................................&
             .................Done"
ENDIF
ENDIF
!-----------------------------------------------------------------------!
!    Perform mesh point leaping	--- 1					!
!-----------------------------------------------------------------------!
WRITE(*,20) "===========================================================&
             ====================="
WRITE(*,20) "Start to perform point leapping............................&
             ....................."
II_toLP=II_PntOvrlp
JJ_toLP=JJ_PntOvrlp
KK_toLP=KK_PntOvrlp
X_toLP=X_Overlap
Y_toLP=Y_Overlap
Z_toLP=Z_Overlap
CALL perform_PntLeap
CALL output_Grid(II_PntLP,JJ_PntLP,KK_PntLP,&
               X_Leap,Y_Leap,Z_Leap,"_PNTLP")
!-----------------------------------------------------------------------!
!   Perform mesh block splitting --- 1	    	 		        !
!-----------------------------------------------------------------------!
WRITE(*,20) "===========================================================&
             ====================="
WRITE(*,20) "Start to perform block splitting...........................&
             ....................."
II_toSplit=II_PntLP
JJ_toSplit=JJ_PntLP
KK_toSplit=KK_PntLP
X_toSplit=X_Leap
Y_toSplit=Y_Leap
Z_toSplit=Z_Leap
CALL perform_BlkSplit
WRITE(*,20) "===========================================================&
             ====================="
CALL write2plt
WRITE(*,20) "===========================================================&
             ====================="
WRITE(*,20) "Start to generate GRIDCON..................................&
             ....................."
CALL gridcon_Gen
WRITE(*,20) "===========================================================&
             ====================="
WRITE(*,20) "Start to generate BCTYPE...................................&
             ....................."
CALL bctype_Gen
WRITE(*,20) "===========================================================&
             ====================="
WRITE(*,20) "Start to generate INIDAMP..................................&
             ....................."
CALL initdamp_Gen
WRITE(*,20) "===========================================================&
             ====================="
WRITE(*,20) "Start to generate run_Case.sh..............................&
             ....................."
CALL run_Gen
STOP
END
