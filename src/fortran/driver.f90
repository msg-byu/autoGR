PROGRAM lat_id_driver
  use find_kgrids, only: grid_selection, find_grids
  use kpointgeneration, only: generateIrredKpointList, mapKptsIntoFirstBZ
  use vector_matrix_utilities, only: matrix_inverse, volume
  use num_types
  ! use fortpy, only: pysave, fpy_read_f, fpy_read
  implicit none

  real(dp) :: lat_vecs(3,3), grid(3,3), offset(3), r_vecs(3,3), point(3), Rinv(3,3)
  real(dp), allocatable :: grids(:,:,:)
  real(dp) :: kpd
  integer :: i, nkpts
  real(dp), pointer :: IRKps(:,:)
  integer, pointer :: weights(:)

  character(300) :: line

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

  call matrix_inverse(transpose(lat_vecs),r_vecs)

  nkpts = int(volume(r_vecs(:,1),r_vecs(:,2),r_vecs(:,3))*kpd)

  call find_grids(lat_vecs,nkpts,grids)
  call grid_selection(lat_vecs, grids, offset, grid)

  call generateIrredKpointList(grid,r_vecs,offset,IRKps,weights)
  call mapKptsIntoFirstBZ(r_vecs,IRKps)
  
  open(4,file="KPOINTS")
  write(4,*) "Our new kpoint method."
  write(4,*) size(IRKps,1)
  write(4,*) "Fractional"
  call matrix_inverse(r_vecs,Rinv)
  do i = 1,size(IRKps,1)
     point = matmul(Rinv,IRKps(i,:))
     write(4,'(F16.14,A1,F16.14,A1,F16.14,A1,I5.1)') point(1), " ",point(2), " ",point(3), " ", weights(i)
  end do

  
  ! print *, "choosen grid"
  ! do i=1,3
  !    print *, grid(:,i)
  ! end do
 
end PROGRAM lat_id_driver
