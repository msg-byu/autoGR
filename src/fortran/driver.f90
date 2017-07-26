PROGRAM lat_id_driver
  use find_kgrids, only: grid_selection, find_grids
  use kpointgeneration, only: generateIrredKpointList
  use vector_matrix_utilities, only: matrix_inverse
  use num_types
  ! use fortpy, only: pysave, fpy_read_f, fpy_read
  implicit none

  real(dp) :: lat_vecs(3,3), grid(3,3), offset(3), r_vecs(3,3)
  real(dp), allocatable :: grids(:,:,:)
  real(dp) :: eps
  integer :: kpd, i
  real(dp), pointer :: IRKps(:,:)
  integer, pointer :: weights(:)

  character(300) :: line

  eps = 1E-4
  open(1,file="POSCAR",status="old")

  read(1,'(a300)') line
  read(1,'(a300)') line
  read(1,*) lat_vecs(:,1)
  read(1,*) lat_vecs(:,2)
  read(1,*) lat_vecs(:,3)
  close(1)

  open(2,file="kpgen",status="old")
  read(2,*) kpd
  read(2,*) offset
  close(2)

  call find_grids(lat_vecs,kpd,grids)
  call grid_selection(grids, grid)

  call matrix_inverse(transpose(lat_vecs),r_vecs)
  call generateIrredKpointList(grid,r_vecs,offset,IRKps,weights,eps_=eps)

  open(4,file="KPOINTS")
  write(4,*) "Our new kpoint method."
  write(4,*) size(IRKps,1)
  write(4,*) "Fractional"
  do i = 1,size(IRKps,1)
     write(4,'(F16.14,A1,F16.14,A1,F16.14,A1,I5.1)') IRKps(i,1), " ",IRKps(i,2), " ",IRKps(i,3), " ", weights(i)
  end do

  
  ! print *, "choosen grid"
  ! do i=1,3
  !    print *, grid(:,i)
  ! end do
 
end PROGRAM lat_id_driver
