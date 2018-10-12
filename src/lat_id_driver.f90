PROGRAM lat_id_driver
  use niggli
  use num_types
  implicit none

  real(dp) :: lat_vecs(3,3), niggli(3,3) Nu(3,3), No(3,3), O(3,3)
  integer :: id, Cu(3,3), Co(3,3), i

  character(300) :: line

  open(1,file="POSCAR",status="old")

  read(1,'(a300)') line
  read(1,'(a300)') line
  do i=1,3
     read(1,*) lat_vecs(:,i)
  end do
  close(1)

  call id_cell(lat_vecs,Nu,Cu,O,No,Co,id)
  call reduce_cell(lat_vecs,niggli,Nu)

  print *, "The cell is niggle cell number: ", id

  open(2,file="Reduced",status="new")
  do i=1,3
     write(2,(' ','3F3.8')) niggli(:,i)
  end do
  
end PROGRAM lat_id_driver
