!!<summary>This module finds the srHNFs for a target volume.
!!Wiley Morgan Started 3/17
!!</summary>
Module srHNF_brute
  use num_types
  use combinatorics, only: binomial
  use numerical_utilities
  use vector_matrix_utilities, only: matrix_inverse
  
  implicit none
  private
  public get_srHNFs,find_pg_gens

CONTAINS

  !!<summary>Finds all the possible diagonals of the HNF matrices of a
  !!given size</summary>
  !!<parameter name="detS" regular="true">Cell size, i.e., determinant
  !!of S matrix.</parameter>
  !!<parameter name="diagonals">All possible diagonals.</parameter>
  SUBROUTINE get_HNF_diagonals(detS, diagonals)
    integer, intent(in) :: detS 
    integer, pointer :: diagonals(:,:) 

    integer i, j, id, quotient, status
    integer :: tempDiag(3,detS*3)
    
    id = 0 ! Number of diagonals found
    do i = 1,detS ! Loop over possible first factors
       if (.not. mod(detS,i)==0) cycle
       quotient = detS/i
       do j = 1,quotient  ! Loop over possible second/third factors
          if (.not. mod(quotient,j)==0) cycle
          id = id + 1
          tempDiag(:,id) = (/i, j, quotient/j /) ! Construct the factor triplet
       enddo
    enddo
    allocate(diagonals(3,id),STAT=status)
    if(status/=0) stop "Allocation failed in get_HNF_diagonals, module deriv..."
    diagonals = tempDiag(:,1:id)
  END SUBROUTINE get_HNF_diagonals
  
  !!<summary>This subroutine does a brute force search to find and return the
  !!HNFs preserve the symmetry of the parent lattice.</summary>
  !!<parameter name="parent" regular="True">The parent lattice.</summary>
  !!<parameter name="gens" regular="True">The generators of the point group.</parameter>
  !!<parameter name="volume" regular="True">Volume factor for the HNF.</parameter>
  !!<parameter name="srHNFS" regular="True">Returned list of symmetry preserving HNFs.</parameter>
  !!<parameter name="eps_" regular="True">epsilon for checking equivalence in floating
  !!point arithmetic.</parameter>
  SUBROUTINE get_srHNFs(parent,gens,volume,srHNFs,eps_)
    real(dp), intent(in) :: gens(:,:,:), parent(3,3)
    integer, intent(in) :: volume
    real(dp), intent(in), optional :: eps_
    integer, allocatable, dimension(:,:,:) :: srHNFs

    integer :: nfound, N, HNF(3,3), current_max
    real(dp) :: eps, B(3,3)
    integer, pointer :: d(:,:) => null()
    integer :: i, j, k, l    ! Loop counters
    integer :: status
    integer, allocatable :: temp_srHNFs(:,:,:)
    integer, allocatable :: temp_storage(:,:,:)

    if(.not. present(eps_)) then
       eps = 1e-10_dp
    else
       eps = eps_
    end if

    call get_HNF_diagonals(volume,d)
    N = size(d,2)
    
    current_max = 10
    allocate(temp_srHNFs(3,3,current_max),STAT=status)
    if(status/=0) stop "Failed to allocate memory in get_srHNFs"
    nfound = 0
    do i = 1,N ! Loop over the permutations of the diagonal elements of the HFNs
       do j = 0,d(2,i)-1  ! Look over possible values of row 2, element 1
          do k = 0,d(3,i)-1  ! Ditto for row 3, element 1
             do l = 0,d(3,i)-1  ! Ditto for row 3, element 2
                HNF(:,:) = reshape((/ d(1,i),      j,     k,        &   
                     0, d(2,i),     l,        &   
                     0,      0, d(3,i)  /), (/3,3/))
                B = matmul(parent,HNF)
                if (symm_pres(B,gens,eps_=eps)) then
                   if (nfound == current_max) then
                      current_max = current_max*2
                      allocate(temp_storage(3,3,size(temp_srHNFs,3)),STAT=status)
                      if(status/=0) stop "Failed to allocate memory in get_srHNFs"
                      temp_storage = temp_srHNFs
                      deallocate(temp_srHNFs)
                      allocate(temp_srHNFs(3,3,current_max),STAT=status)
                      if(status/=0) stop "Failed to allocate memory in get_srHNFs"
                      temp_srHNFs = 0
                      temp_srHNFs(:,:,1:nfound) = temp_storage(:,:,1:nfound)
                      deallocate(temp_storage)
                   end if
                   nfound = nfound+1 ! Count the srHNFs
                   temp_srHNFs(:,:,nfound) = HNF
                end if
             enddo;enddo;enddo  ! End loops over values for off-diagonal elements
    enddo ! End loop over all unique triplets of target determinant (volume)
    allocate(srHNFs(3,3,nfound))
    srHNFs(:,:,1:nfound) = temp_srHNFs(:,:,1:nfound)
    if (associated(d)) then
       deallocate(d)
    end if
    if (allocated(temp_storage)) then
       deallocate(temp_storage)
    end if
    deallocate(temp_srHNFs)
  end SUBROUTINE get_srHNFs

  !!<summary>This function determines if a lattice has the desired
  !!symmetries.</summary>
  !!<parameter name="lat" regular="True">The super lattice vectors.</parameter>
  !!<parameter name="pg" regular="True">The poing group operators we want to check for.</parameter>
  !!<parameter name="eps_" regular="True">epsilon for checking equivalence in floating
  !!point arithmetic.</parameter>
  logical function symm_pres(lat,pg,eps_)
    real(dp), intent(in) :: lat(3,3)
    real(dp), intent(in) :: pg(:,:,:)
    real(dp), intent(in), optional :: eps_

    real(dp) :: eps
    integer :: i
    real(dp) :: RB(3,3), N(3,3), latinv(3,3)
    
    if(.not. present(eps_)) then
       eps = 1e-10_dp
    else
       eps = eps_
    end if

    symm_pres = .True.
    do i=1,size(pg,3)
       RB = matmul(pg(:,:,i),lat)
       call matrix_inverse(lat,latinv)
       N = matmul(latinv,RB)
       if (.not. equal(N,nint(N),eps)) then
          symm_pres = .False.
          exit
       end if
    end do
    
  end function symm_pres

  !!<summary>Finds the generators of a given point group.</summary>
  !!<parameter name="pg" regular="True">The point group.</parameter>
  !!<parameter name="gens" regular="True">The generators of the point group.</parameter>
  !!<parameter name="eps" regular="True">epsilon for checking equivalence in floating
  !!point arithmetic.</parameter>
  subroutine find_pg_gens(pg,gens,eps)
    real(dp), pointer, intent(in) :: pg(:,:,:)
    real(dp), allocatable, intent(out) :: gens(:,:,:)
    real(dp), intent(in) :: eps
    
    real(dp), allocatable :: test_gens(:,:,:)
    integer :: n_gens, status, i, j, gens_found
    logical :: found_gens
    integer, allocatable :: combs(:,:)

    gens_found = 0
    n_gens = 2
    found_gens = .False.

    do while ((found_gens .eqv. .False.) .and. (n_gens < 4))
       call combinations(size(pg,3),n_gens,combs)
       allocate(test_gens(3,3,n_gens),STAT=status)
       if(status/=0) stop "Failed to allocate memory in find_pg_gens"
       do i=1,size(combs,2)
          do j=1,n_gens
             test_gens(:,:,j) = pg(:,:,combs(j,i))
          end do
          if (forms_group(test_gens,pg,eps) .eqv. .True.) then
             allocate(gens(3,3,n_gens),STAT=status) 
             if(status/=0) stop "Failed to allocate memory in find_pg_gens"
             gens(:,:,:) = test_gens(:,:,:)
             found_gens = .True.
             exit
          end if
       end do

       n_gens = n_gens + 1
       deallocate(test_gens)
    end do

    if (allocated(combs)) then
       deallocate(combs)
    end if
  end subroutine find_pg_gens

  !!<summary>Determines if a set of generators forms the group.</summary>
  !!<parameter name="gens" regular="True">The test set of generators.</parameter>
  !!<parameter name="pg" regular="True">The pg we want to be able to generate.</parameter>
  !!<parameter name="eps" regular="True">epsilon for checking equivalence in floating
  !!point arithmetic.</parameter>
  logical function forms_group(gens,pg,eps)
    real(dp), intent(in) :: gens(:,:,:), pg(:,:,:)
    real(dp), intent(in) :: eps
    
    integer :: elements, i, j, k, nf, temp_int, status
    real(dp), allocatable :: test_set(:,:,:)
    real(dp) :: test(3,3)
    logical :: growing, in_group


    forms_group = .True.
    allocate(test_set(3,3,size(pg,3)),STAT=status)
    if(status/=0) stop "Failed to allocate memory in forms_group"
    test_set = 0
    elements = 0

    do i=1,size(gens,3)
       test_set(:,:,i) = gens(:,:,i)
       elements = elements + 1
    end do

    do i=1, size(gens,3)
       do j=1,size(gens,3)
          test = matmul(gens(:,:,i),gens(:,:,j))
          in_group = .False.
          do k=1,size(test_set)
             if (equal(test,test_set(:,:,k),eps)) then
                in_group = .True.
                exit
             end if
          end do
          if (in_group .eqv. .False.) then
             elements = elements + 1
             if (elements > size(test_set,3)) then
                forms_group = .False.
                exit
             end if
             test_set(:,:,elements) = test
          end if
       end do
       if (forms_group .eqv. .False.) then
          exit
       end if
    end do

    growing = .True.
    do while (growing .eqv. .True.)
       nf = 0
       temp_int = elements
       do i=1, size(gens,3)
          do j=1,temp_int
             test = matmul(gens(:,:,i),test_set(:,:,j))
             in_group = .False.
             do k=1,size(test_set)
                if (equal(test,test_set(:,:,k),eps)) then
                   in_group = .True.
                   exit
                end if
             end do
             if (in_group .eqv. .False.) then
                elements = elements + 1
                nf = nf + 1
                if (elements > size(test_set,3)) then
                   forms_group = .False.
                   exit
                end if
                test_set(:,:,elements) = test
             end if
          end do
          if (forms_group .eqv. .False.) then
             exit
          end if
       end do

       if ((nf == 0) .or. (forms_group .eqv. .False.)) then
          growing = .False.
       end if
    end do

    if (size(test_set,3) /= size(pg,3)) then
       forms_group = .False.
    else
       do i=1,size(pg,3)
          in_group = .False.
          do j=1,size(test_set)
             if (equal(pg(:,:,i),test_set(:,:,j),eps)) then
                forms_group = .True.
                exit
             else
                forms_group = .False.
             end if
          end do
          if (forms_group .eqv. .False.) then
             exit
          end if
       end do
    end if
    deallocate(test_set)
  end function forms_group
  
  !!<summary>Finds all pair or triplet combinations of numbers 1 to n.</summary>
  !!<parameter name="n" regular="True">The largest number in the pairs.</parameter>
  !!<parameter name="comb" regular="True">The resulting combinations.</parameter>
  !!<parameter name="r" regular="True">2 if for pairs, 3 if for triplets.</parameter>
  subroutine combinations(n,r,comb)
    integer, intent(in) :: n, r
    integer, allocatable, intent(out) :: comb(:,:)

    integer :: n_combs, status, i, j, k, combs_found

    if (allocated(comb)) then
       deallocate(comb)
    end if
    if (r == 2) then
    
       n_combs = binomial(n,2)
       combs_found = 0
       allocate(comb(3,n_combs),STAT=status)
       if(status/=0) stop "Failed to allocate memory in combinations"
       do i = 1,n-1
          do j = i+1, n
             combs_found = combs_found + 1
             comb(:,combs_found) = (/i,j,0/)
          end do
       end do
    elseif (r==3) then
       n_combs = binomial(n,3)
       combs_found = 0
       allocate(comb(3,n_combs),STAT=status)
       if(status/=0) stop "Failed to allocate memory in combinations"
       do i = 1,n-2
          do j = i+1, n-1
             do k = j+1, n
                combs_found = combs_found + 1
                comb(:,combs_found) = (/i,j,k/)
             end do
          end do
       end do
    else
       stop "combinations only works for pair and triplets (r=2 or 3)"
    end if
  end subroutine combinations
end Module srHNF_brute
