!!<summary>Generates the symmetry preserving HNFs for the different
!!crystal lattices.</summary>
Module sp_hnfs
  implicit none
  private
  public sc

CONTAINS

  !!<summary>Finds the symmetry preserving HNFs for the simple cubic
  !!lattice with determinant n. Assuming the basis of A =
  !![[1,0,0],[0,1,0],[0,0,1]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  SUBROUTINE sc(n,spHNFs)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs
    integer :: temp_HNFs(3,3,1)
    
    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0

    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((a==c) .and. (c==f)) then
          b = 0
          d = 0
          e = 0
          nhnfs = nhnfs + 1
          temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((a==c) .and. (real(f)/real(a)==2)) then
          b = 0
          d = a
          e = a
          nhnfs = nhnfs + 1          
          temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((c==f) .and. (real(f)/real(a)==2)) then
          b = a
          d = a
          e = 0
          nhnfs = nhnfs + 1          
          temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       end if
    end do
    allocate(spHNFs(3,3,nhnfs))

    spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    
  end SUBROUTINE sc

  !!<summary>Finds the symmetry preserving HNFs for the face centered
  !!cubic lattice with determinant n. Assuming the basis of A =
  !![[0,1,1],[1,0,1],[1,1,0]].</summary>  
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  SUBROUTINE fcc(n,spHNFs)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs
    integer :: temp_HNFs(3,3,1)
    
    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0

    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((a==c) .and. (a==f)) then
          b = 0
          d = 0
          e = 0
          nhnfs = nhnfs + 1          
          temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((c==f) .and. ((real(f)/real(a) == 2) .or. (real(f)/real(a) == 4))) then
          b = a
          d = a
          e = 0
          nhnfs = nhnfs + 1          
          temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       end if
    end do
    allocate(spHNFs(3,3,nhnfs))

    spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    
  end SUBROUTINE fcc

  !!<summary>Finds the symmetry preserving HNFs for the body centered
  !!cubic lattice with determinant n. Assuming the basis of A =
  !![[-1,1,1],[1,-1,1],[1,1,-1]].</summary>  
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  SUBROUTINE bcc(n,spHNFs)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs
    integer :: temp_HNFs(3,3,1)
    
    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0

    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((a==c) .and. (a==f)) then
          b = 0
          d = 0
          e = 0
          nhnfs = nhnfs + 1          
          temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((a==c) .and. (real(f)/real(a)==2)) then
          d = 0
          d = a
          e = a
          nhnfs = nhnfs + 1          
          temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((a==c) .and. (real(f)/real(a) ==4)) then
          b = 0
          d = 3*a
          e = 3*a
          nhnfs = nhnfs + 1          
          temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       end if
    end do
    allocate(spHNFs(3,3,nhnfs))

    spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    
  end SUBROUTINE bcc

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
  






end Module sp_hnfs
