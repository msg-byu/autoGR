PROGRAM driver
  use find_kgrids, only: grid_selection, find_grids
  use kpointgeneration, only: generateIrredKpointList, mapKptsIntoFirstBZ
  use symmetry, only: bring_into_cell
  use vector_matrix_utilities, only: matrix_inverse
  use num_types
  ! use fortpy, only: pysave, fpy_read_f, fpy_read
  implicit none

  real(dp) :: lat_vecs(3,3), grid(3,3), offset(3), r_vecs(3,3), point(3), Rinv(3,3)
  real(dp), allocatable :: grids(:,:,:), offsets(:,:)
  real(dp) :: eps
  integer :: kpd, i
  real(dp), pointer :: IRKps(:,:)
  integer, pointer :: weights(:)

  character(300) :: line

  eps = 1E-6
  open(1,file="POSCAR",status="old")

  read(1,'(a300)') line
  read(1,'(a300)') line
  read(1,*) lat_vecs(:,1)
  read(1,*) lat_vecs(:,2)
  read(1,*) lat_vecs(:,3)
  close(1)

  open(2,file="kpgen",status="old")
  read(2,*) kpd
  ! read(2,*) offset
  close(2)

  call matrix_inverse(transpose(lat_vecs),r_vecs)
  
  call find_grids(lat_vecs,kpd,grids,offsets,eps_=eps)
  print *, "grids to search", size(grids,3)
  call grid_selection(grids, grid, offset, r_vecs, offsets,eps_=eps)
  
  call generateIrredKpointList(grid,r_vecs,offset,IRKps,weights,eps_=eps)
  ! call mapKptsIntoFirstBZ(r_vecs,IRKps,eps_=eps)

  open(4,file="KPOINTS")
  write(4,'(A29,I6.1,A14)') "Our new kpoint method. Found ", sum(weights), " total points."
  write(4,*) size(IRKps,1)
  write(4,*) "Fractional"
  call matrix_inverse(r_vecs,Rinv)
  do i = 1,size(IRKps,1)
     point = IRKps(i,:)
     call bring_into_cell(point,Rinv,r_vecs,eps)
     point = matmul(Rinv,IRKps(i,:))
     write(4,'(F16.14,A1,F16.14,A1,F16.14,A1,F5.1,A3,I5)') point(1), " ",point(2), " ",point(3), " ", real(weights(i)), " ! ", i
  end do

  
  ! print *, "choosen grid"
  ! do i=1,3
  !    print *, grid(:,i)
  ! end do
 
end PROGRAM driver
