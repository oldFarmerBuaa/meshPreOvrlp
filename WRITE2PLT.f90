subroutine write2plt

use system_var

implicit none

integer imb,mb,imax,jmax,ij,il,jl,ijmax,ij_start,ij_end
character(len=14) fn
character(len=10) fn_prefix
character(len=18) fn_plt
real(8),dimension(:),allocatable :: x,y

!write(*,*) "ls GRIDBLO*|wc -l"
!call system("ls GRIDBLO*|wc -l")
!call system("ls GRIDBLO*.plt|wc -l")
!write(*,*) "Please input number of total blocks(Take smaller number):"
mb=n_BlkFinal
fn_prefix='GRIDBLO_BC'
write(*,*) 'mb=',mb,'filename=',fn_prefix
do imb=1,mb
  write(*,*) 'i_Block=:',imb
  write(fn,'(A,I4.4)')fn_prefix,imb
  open(10,file=fn,status='old')
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*) jmax,imax
  ijmax=imax*jmax
  allocate (x(ijmax),y(ijmax))
  ij_start=1
  ij_end=ijmax
  read(10,*) (x(ij),ij=ij_start,ij_end)
  read(10,*) (y(ij),ij=ij_start,ij_end)
  fn=trim(fn)
  write(fn_plt,'(A,A)')fn,'.plt'
  open(20,file=fn_plt,status='unknown')
  write(20,*) 'TITLE = "',fn,'"'
  write(20,*) 'VARIABLES = ','"X"',',','"Y"'
  write(20,'(A5,A2,I6,A3,I6,A20)')'ZONE ','I=',imax,',J=',jmax,',DATAPACKING=POINT'
  do ij=ij_start,ij_end
    write(20,*) x(ij),y(ij)
  enddo
  close(10)
  close(20)
  deallocate(x,y)
enddo
write(*,*) 'mb=',mb,'filename=',fn_prefix
write(*,*) 'Done.'
end
