!-----------------------------------------------------------------------!
! Variable name rule: <main property>_<other properties>
!			e.g., min_LvOffUniformTot,n_Point_LP_Loc
! Variable prefix:
!		i_ : index of sth.
!		n_ : total number of sth.
!		II_,JJ_,KK_ : structured mesh max index I,J,K
!		Lv_: grid evel of mesh block, smaller = finer mesh
!		tmp_ : temperal values for Int,Real,Char
! Variable suffix:
!		*Pre : of original mesh
!		*Loc : in local region
!		*Glb : in global region
!
!-----------------------------------------------------------------------!
! Reading GRIDCON_Pre
MODULE READ_GRIDCON_Pre
IMPLICIT NONE
CHARACTER(LEN=2) cas_Dim    ! case dimension 2D or 3D
CHARACTER(LEN=6) mesh_Type  ! mesh type to read: cfx4,fluent,cgns...
INTEGER n_BlkSupPre,n_BlkGlbPre
INTEGER max_BlkGlbPre
INTEGER min_LvOffUniformTot
INTEGER n_PointOverlap
INTEGER,DIMENSION(:),ALLOCATABLE :: n_BlkPre,max_GlueInBlkPre
INTEGER,DIMENSION(:),ALLOCATABLE :: II_Blk,JJ_Blk,KK_Blk 
INTEGER,DIMENSION(:,:),ALLOCATABLE :: Lv_OffUniformTot,Lv_Uniform
INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: isidtp_Pre,&
                                        Lv_OffSet,Lv_Coarsen,&
                                        Lv_OffUniform
! cas_Dim    : case dimension 2D or 3D
! mesh_Type  : mesh type to read: cfx4,fluent,cgns...
! min_LvOffUniformTot: min final grid level
! n_PointOverlap : total overlaping points
! n_BlkSupPre: total superblocks
! n_BlkGlbPre: total blocks in all superblocks
! n_BlkPre(i_BlkSupPre) : total blocks in each superblock
! max_BlkGlbPre: max number of blocks among all uperblock
! max_GlueInBlkPre(i_BlkSupPre) : max glue in blocks of each superblock
! II/JJ/KK_Blk : Max idx of block matrix in each superblock
! Lv_OffUniformTot : total level of three directions
! isidtp_Pre : boundary conditions
! Lv_OffSet  : grid level of origianl block compared to uniform mesh
! Lv_Coarsen : level to be coarsened
! Lv_OffUniform : final level compared to uniform mesh
END MODULE

!-----------------------------------------------------------------------!
! Reading original mesh
MODULE read_mesh_original
IMPLICIT NONE
INTEGER n_PathchPre,n_GluePre,n_ElmPre
INTEGER max_PointPre,max_PointLocPre
INTEGER,DIMENSION(:),ALLOCATABLE :: n_PointPre
INTEGER,DIMENSION(:,:),ALLOCATABLE :: II_Pre,JJ_Pre,KK_Pre,&
                                      n_PointLocPre
REAL(8),DIMENSION(:,:),ALLOCATABLE :: X_Pre,Y_Pre,Z_Pre
! n_PathchPre : number of patches in each superblock
! n_GluePre   : numnber of glues in each superblock
! n_ElmPre    : number of elements in each superblock
! max_PointPre : max number of points among all superblocks
! max_PointLocPre: max point number among all blocks
! n_PointPre(i_BlkSupPre) : total points in each superblock
! II/JJ/KK_Pre : Index of I,J,K in each block
! n_PointLocPre : number of point in each block
! X/Y/Z_Pre: Initial X/Y/Z read from mesh file
END MODULE

!-----------------------------------------------------------------------!
! System variables: I/O, tmp, memory ratio, logs
MODULE system_var
IMPLICIT NONE
INTEGER tmp_Int,tmp_Int1,tmp_Int2,tmp_Int3,CFL
INTEGER ratio_PointAftPre,np,flag_Mltstp
INTEGER n_BlkFinal,Lv_Mltstp,flag_CalTmStpIdx,n_PntFnlTot,max_Proc
REAL(8) tmp_Real,tmp_Real1,tmp_Real2,tmp_Real3
CHARACTER(LEN=20) tmp_Char
CHARACTER(LEN=18) LOGFILE
CHARACTER(LEN=80) message
CHARACTER(LEN=30) FILENAME
INTEGER,DIMENSION(:),ALLOCATABLE :: II_Final,JJ_Final,KK_Final,n_PntFinal
INTEGER,DIMENSION(:),ALLOCATABLE :: idx_Proc,idx_Mltstp
INTEGER,DIMENSION(:,:),ALLOCATABLE :: Lv_OffUnfrmFnl,isidtp_Final
END MODULE
! flag_Mltstp : flag of multi-time-step or not, flag_Mltstp=1: use multi-time-step
! n_BlkFinal: Final No. blocks after all processing.
! flag_CalTmStpIdx: flag of calulate time step index of GRIDCON, 1 = yes, 0 = no
! n_PntFnlTot : Total final points of all blocks
! max_Proc : Max processor index assigned.
! Lv_OffUnfrmFnl: (3,n_BlkFinal), Lv_OffUniform of final blocks.
! isidtp_Final: (6,n_BlkFinal), boundary conditions of final blocks.
!-----------------------------------------------------------------------!
! Block point leaping
MODULE point_leaping
IMPLICIT NONE
INTEGER Lv_CoarsenGlb
INTEGER,DIMENSION(:),ALLOCATABLE :: n_PntAftLp
INTEGER,DIMENSION(:,:),ALLOCATABLE :: n_PointLP
INTEGER,DIMENSION(:,:),ALLOCATABLE :: II_toLP,JJ_toLP,KK_toLP
INTEGER,DIMENSION(:,:),ALLOCATABLE :: II_PntLP,JJ_PntLP,KK_PntLP
REAL(8),DIMENSION(:,:),ALLOCATABLE :: X_toLP,Y_toLP,Z_toLP
REAL(8),DIMENSION(:,:),ALLOCATABLE :: X_Leap,Y_Leap,Z_Leap
! Lv_CoarsenGlb: Level of mesh global coarsen
! n_PntAftLp : Total points after leaping in each super block.
! II/JJ/KK_toLP: index of I,J,K before leaping in each block
! II/JJ/KK_PntLP: index of I,J,K after leaping in each block
! X_toLP,Y_toLP,Z_toLP: coordinates before leaping
! X_Leap,Y_Leap,Z_Leap: coordinates after leaping
!
END MODULE

!-----------------------------------------------------------------------!
! Block overlap
MODULE block_overlaping
IMPLICIT NONE
INTEGER f_Ovrlp1Or2Dir
INTEGER,DIMENSION(:),ALLOCATABLE   :: n_PntToOvrlp,n_PntAftOvrlp
INTEGER,DIMENSION(:,:),ALLOCATABLE   :: n_PntAftOvrlpBfBlk
INTEGER,DIMENSION(:,:),ALLOCATABLE   :: n_PntInBlk,n_PntBfBlk
INTEGER,DIMENSION(:,:),ALLOCATABLE :: II_toOvrlp,JJ_toOvrlp,KK_toOvrlp
INTEGER,DIMENSION(:,:),ALLOCATABLE :: II_PntOvrlp,JJ_PntOvrlp,KK_PntOvrlp
INTEGER,DIMENSION(:,:),ALLOCATABLE :: n_PntOvrlpLoc
INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: dir_Swipe
REAL(8),DIMENSION(:,:),ALLOCATABLE :: X_toOvrlp,Y_toOvrlp,Z_toOvrlp
REAL(8),DIMENSION(:,:),ALLOCATABLE :: X_Overlap,Y_Overlap,Z_Overlap
! f_Ovrlp1Or2Dir : Flag of overlap pattern, 1=One-direction, 2=Two-directions
! n_PntToOvrlp : No. points before overlapping in each super block.
! n_PntAftOvrlp : No. points after overlapping in each super block.
! n_PntAftOvrlpBfBlk(J,K): No. pnts after overlap  before block J.
! n_PntInBlk: No. points in block, used as interval points number array.
! II/JJ/KK_PntLP: index of I,J,K after overlaping in each block
! n_PntOvrlpLoc: No. points in each block after overlapping.
! dir_Swipe : Direction of swipe in overlapping in each block (6,n_Blk,n_BlkSup)
END MODULE

!-----------------------------------------------------------------------!
! Block splitting
MODULE block_splitting
IMPLICIT NONE
INTEGER max_N_Step,flag_MoreSplit,n_MoreSplit
INTEGER,DIMENSION(:,:),ALLOCATABLE :: n_PntEqvLoc,n_toSplitLoc
INTEGER,DIMENSION(:,:),ALLOCATABLE :: n_Split
INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: n_step,Lv_DiffVsMin
INTEGER,DIMENSION(:,:),ALLOCATABLE :: II_toSplit,JJ_toSplit,KK_toSplit
REAL(8),DIMENSION(:),ALLOCATABLE :: min_dxyz,min_dt
REAL(8),DIMENSION(:,:,:),ALLOCATABLE :: dt
REAL(8),DIMENSION(:,:),ALLOCATABLE :: X_toSplit,Y_toSplit,Z_toSplit
INTEGER,PARAMETER :: n_PntOvrlpInSplit=6
END MODULE
! max_N_Step : max No. steps needed to march to Time = 1 among all n_step
! flag_MoreSplit: flag of more split than initial balance or not, 1 = more split, 0 = no more
! n_MoreSplit: No. of extra splits activated by flag_MoreSplit=1.
! n_PntEqvLoc: Equivalent points in each block n_PntEqvLoc(J,K)
! n_toSplitLoc: local No. points in block (J,K)
! n_Split    : No. split in each block n_Split(J,K)
! n_step     : No. steps in each block needed to march to Time = 1, n_step(1:3,J,K)
! Lv_DiffVsMin: Lv_Diff vs min Lv 
! min_dxyz   : min dxyz in each super block min_dxyz(K)
! min_dt     : min dt in each super block min_dt(K)
! dt         : dt of each block dt(1:3,J,K)
