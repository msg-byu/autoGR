PROGRAM lat_id_driver
    use sp_hnfs
    use rational_mathematics, only: SmithNormalForm
    use num_types

    implicit none

    integer :: n, at(1), n_irr, nhnfs, i, err, j,k
    integer :: r_low, r_high
    CHARACTER(len=32) :: arg1,arg2
    integer, allocatable :: HNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer :: Co(3,3), Cu(3,3), S(3,3), L(3,3), R(3,3)
    real(dp) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(1,3)
    real(dp) :: grid(3,3), best_offset(3)
    logical :: all_hnfs
    !Get command line args
    call getarg(1,arg1)
    call getarg(2,arg2)
    read(arg1,*) r_low
    read(arg2,*) r_high

    ! This is just some basic setup stuff which you should never chance
    ! since it's not effecting the code at all but it has to be here for
    ! the code to run.
    Co = 0; Cu = 0; No = 0.0_dp; Nu = 0.0_dp; O = 0.0_dp; U = 0.0_dp; offsets = 0.0_dp;
    at = 1
    all_hnfs = .True.
    allocate(B_vecs(3,1))
    B_vecs = 0.0_dp
    do i=1,3
    Co(i,i) = 1; Cu(i,i) = 1
    No(i,i) = 1.0_dp; O(i,i) = 1.0_dp; Nu(i,i) = 1.0_dp; U(i,i) = 1.0_dp
    end do
    k=0
    ! Change the upper and lower limits on this to whatever you want.
    do n = r_low,r_high
    print *, "Running HNFs of determinant: ", n
    ! Change this function call to whichever case from sp_hnfs you want.
    if(Mod(n,1000) == 0) then
        print *, "Running HNFs of determinant: ", n
    endif
    call hex_12(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, HNFs, grid, &
        n_irr, nhnfs, all_hnfs_=all_hnfs)
    print *, "nhnfs", nhnfs
    if (nhnfs > 0) then
        k = k+nhnfs
        err = 0
        do i =1, nhnfs
        call SmithNormalForm(HNFs(:,:,i),L, S, R, err_=err)
        if (err == 1) then
            print *, "overflow ",n
            print *, "HNF"
            do j=1,3
            print *, HNFs(j,:,i)
            end do
        end if
        end do
    end if
    end do 
print *, "total number of HNFS", k

end PROGRAM lat_id_driver
