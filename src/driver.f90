PROGRAM lat_id_driver
  use find_kgrids, only: find_grid
  use kpointgeneration, only: generateIrredKpointList, mapKptsIntoBZ
  use control_file, only: get_inputs
  use vector_matrix_utilities, only: matrix_inverse, determinant, minkowski_reduce_basis
  use numerical_utilities, only: equal
  use num_types

  implicit none

  real(dp) :: lat_vecs(3,3), grid(3,3), offset(3), r_vecs(3,3), reduced_R(3,3)
  real(dp) :: Rinv(3,3), point(3), eps, best_offset(3)
  logical :: find_offset, min_kpts
  integer :: nkpts, i, symm_flag
  integer, allocatable :: at(:)
  real(dp), pointer :: IRKps(:,:)
  real(dp), allocatable :: B_vecs(:,:)
  integer, pointer :: weights(:)

  call get_inputs(nkpts, lat_vecs, at, B_vecs, offset, find_offset, symm_flag, &
       min_kpts, eps)
  call matrix_inverse(transpose(lat_vecs),r_vecs)
  call minkowski_reduce_basis(r_vecs, reduced_R, eps)
  call find_grid(lat_vecs, nkpts, B_vecs, at, offset, find_offset, grid, best_offset, &
       symm_flag_=symm_flag, min_kpts_=min_kpts, eps_=eps)

  call generateIrredKpointList(lat_vecs, B_vecs, at, grid, reduced_R, best_offset, &
       IRKps, weights, reps_=eps, symm_=symm_flag)

  call mapKptsIntoBZ(r_vecs, IRKps, eps)

  open(4,file="KPOINTS")
  write(4,'(A69, I12, A1, I6)')"Dynamic K-point generation: number of irreducible / number reducible: ", size(IRKps,1), '/', sum(weights)
  write(4,*) size(IRKps,1)
  write(4,*) "Fractional"

  call matrix_inverse(r_vecs,Rinv)
  do i = 1,size(IRKps,1)
     point = matmul(Rinv,IRKps(i,:))
     write(4,'(F32.14,A1,F32.14,A1,F32.14,A1,I5.1)') point(1), " ",point(2), " ",point(3), " ", weights(i)
  end do

end PROGRAM lat_id_driver
