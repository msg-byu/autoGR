PROGRAM lat_id_driver
  use find_kgrids, only: find_grid
  use kpointgeneration, only: generateIrredKpointList, mapKptsIntoBZ
  use control_file, only: get_inputs
  use vector_matrix_utilities, only: matrix_inverse, determinant, minkowski_reduce_basis
  use numerical_utilities, only: equal
  use num_types
  ! use fortpy, only: pysave, fpy_read_f, fpy_read
  implicit none

  real(dp) :: lat_vecs(3,3), grid(3,3), offset(3), r_vecs(3,3), reduced_R(3,3)
  real(dp) :: Rinv(3,3), point(3), reps, aeps, best_offset(3), det
  logical :: find_offset
  integer :: nkpts, i
  integer, allocatable :: at(:)
  real(dp), pointer :: IRKps(:,:), B_vecs(:,:)
  integer, pointer :: weights(:)

  call get_inputs(nkpts, lat_vecs, at, B_vecs, offset, find_offset, reps, aeps)
  call matrix_inverse(transpose(lat_vecs),r_vecs)
  call find_grid(lat_vecs, nkpts, B_vecs, at, offset, find_offset, grid, best_offset, &
       eps_=reps)

  det = determinant(grid)
  ! if (((det<eps) .and. (det >1E-10)) .or. (equal(det,eps,eps))) then
  !    call generateIrredKpointList(lat_vecs, B_vecs, at, grid, r_vecs, best_offset, &
  !         IRKps, weights, eps_=(eps**2))
  ! else
  ! print *, "offset", best_offset
  ! print *, "grid"
  ! do i=1,3
  !    print *, grid(i,:)
  ! end do
  call generateIrredKpointList(lat_vecs, B_vecs, at, grid, r_vecs, best_offset, &
       IRKps, weights, reps)
  ! end if
  call mapKptsIntoBZ(r_vecs, IRKps, reps)

  open(4,file="KPOINTS")
  ! write(4,'("Our new kpoint method. ",i6)') sum(weights)
  write(4,'(A69, I12, A1, I6)')"Dynamic K-point generation number of irreducible / number reducible: ", size(IRKps,1), '/', sum(weights)
  write(4,*) size(IRKps,1)
  write(4,*) "Fractional"
  ! write(4,'(A4)')"cart"
  call matrix_inverse(r_vecs,Rinv)
  do i = 1,size(IRKps,1)
     point = matmul(Rinv,IRKps(i,:))
     write(4,'(F16.14,A1,F16.14,A1,F16.14,A1,I5.1)') point(1), " ",point(2), " ",point(3), " ", weights(i)
  end do

end PROGRAM lat_id_driver
