PROGRAM lat_id_driver
  use lat_types
  use num_types
  use fortpy, only: pysave, fpy_read_f, fpy_read
  implicit none

  real(dp) :: u_vecs(3,3), b_vecs(3,3)
  real(dp) :: eps
  integer :: lat_id

  character(300) :: line

  eps = 1E-4
  ! open(1,file="POSCAR",status="old")

  ! read(1,'(a300)') line
  ! read(1,'(a300)') line
  ! read(1,*) u_vecs(:,1)
  ! read(1,*) u_vecs(:,2)
  ! read(1,*) u_vecs(:,3)
  ! close(1)

  call fpy_read_f("identify_lattice_vecs.in.tric2","#",u_vecs)
  call fpy_read("identify_lattice_lat_id.out.tric2","#",lat_id)
  print *, "lat_id",lat_id
  call canonical_basis(lat_id,u_vecs,b_vecs)
  call pysave(b_vecs,"canonical_basis.out.tric2")

  ! print *,"lat_id", lat_id
  
end PROGRAM lat_id_driver
