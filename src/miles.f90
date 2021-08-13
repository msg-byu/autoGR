PROGRAM driver

  use sp_hnfs, only: so_32
  use num_types
  
  implicit none

  integer :: n, i, j
  integer, allocatable :: spHNFs(:,:,:)
  real(dp), pointer :: B_vecs(:,:)
  integer :: at(1)
  integer :: Co(3,3), Cu(3,3)
  real(dp) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(1,3)
  integer :: n_irr, nhnfs
  real(dp) :: grid(3,3), best_offset(3)
  logical :: all_hnfs

  n = 100
  allocate(B_vecs(1,3))
  B_vecs(1,:) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  offsets(1,:) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  at(1) = 0
  all_hnfs = .True.

  Co(1,:) = (/ 1, 0, 0 /)
  Co(2,:) = (/ 0, 1, 0 /)
  Co(3,:) = (/ 0, 0, 1 /)

  Cu(1,:) = (/ 1, 0, 0 /)
  Cu(2,:) = (/ 0, 1, 0 /)
  Cu(3,:) = (/ 0, 0, 1 /)
  
  No(1,:) = (/1.0_dp, 0.0_dp, 0.0_dp /)
  No(2,:) = (/0.0_dp, 1.0_dp, 0.0_dp /)
  No(3,:) = (/0.0_dp, 0.0_dp, 1.0_dp /)

  Nu(1,:) = (/1.0_dp, 0.0_dp, 0.0_dp /)
  Nu(2,:) = (/0.0_dp, 1.0_dp, 0.0_dp /)
  Nu(3,:) = (/0.0_dp, 0.0_dp, 1.0_dp /)

  O(1,:) = (/1.0_dp, 0.0_dp, 0.0_dp /)
  O(2,:) = (/0.0_dp, 1.0_dp, 0.0_dp /)
  O(3,:) = (/0.0_dp, 0.0_dp, 1.0_dp /)
  
  U(1,:) = (/1.0_dp, 0.0_dp, 0.0_dp /)
  U(2,:) = (/0.0_dp, 1.0_dp, 0.0_dp /)
  U(3,:) = (/0.0_dp, 0.0_dp, 1.0_dp /)

  call so_32(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, all_hnfs_=all_hnfs)

  open(3, file="HNFs.out")
  do i= 1,nhnfs
     write(3,'(6i6)') spHNFs(1,1,i), spHNFs(2,1,i), spHNFs(2,2,i), spHNFs(3,1,i), spHNFs(3,2,i), spHNFs(3,3,i)
  end do
  close(3)
     
end PROGRAM driver
