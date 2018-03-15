PROGRAM driver
  use srHNF_brute
  use symmetry, only: get_lattice_pointGroup
  use vector_matrix_utilities
  use num_types
  use utilities
  use io_utils
  use numerical_utilities

  implicit none

  integer :: i, sfound, nfound, status, j
  integer, allocatable :: srHNFs(:,:,:), isrHNFs(:,:,:), temp_storage(:,:,:)
  real(dp), pointer :: pg(:,:,:)
  real(dp), allocatable :: gens(:,:,:)
  
  integer nMin, nMax ! Numbers of various things
  integer k ! Number of colors/label types (i.e., binary, ternary, etc.)
  integer LatDim ! 2D or 3D parent lattice?
  integer nD ! Number of sites in the basis (i.e., number of points in the multilattice)
  real(dp), pointer :: d(:,:) => null()
  character(1) :: latTyp
  real(dp)  eps, tstart, tend
  real(dp) :: parLV(3,3), invLV(3,3), abgens(3,3)

  character(80) title, fname
  logical fullLab,concCheck
  integer, pointer :: cRange(:,:)
  integer, pointer :: label(:,:)
  integer, pointer :: digit(:)
  integer, pointer :: equivalencies(:)

  if (iargc()>=1) then
     call getarg(1,fname)
  else
     fname = "struct_enum.in"
  endif
  call read_input(title,LatDim,parLV,nD,d,k,equivalencies,nMin,nMax,eps&
       &,fullLab,label,digit,fname,cRange,concCheck) ! Read in parent lattice vectors, etc.
  
  ! lat = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
  ! lat = reshape((/0.5,0.5,0.0,0.5,0.0,0.5,0.0,0.5,0.5/),(/3,3/))
  ! lat = reshape((/0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5,0.5/),(/3,3/))
  ! vol = 1000

  call cpu_time(tstart)
  ! Need to replace get_lattice_pointGroup with get space group
  call get_lattice_pointGroup(parLV,pg)
  call find_pg_gens(pg,gens,1e-10_dp)


  ! print *, "pg size", size(pg,3)
  ! print *, "lat vecs", parLV

  do i=1, size(gens,3)
     ! print *, "cbas gen", i
     ! do j=1,3
     !    print *, gens(j,:,i)
     ! end do
     print *, "latbas gen", i
     call matrix_inverse(parLV,invLV)
     abgens = matmul(invLV,matmul(gens(:,:,i),parLV))
     do j=1,3
        print *, abgens(j,:)
     end do
  end do
  
end PROGRAM driver
