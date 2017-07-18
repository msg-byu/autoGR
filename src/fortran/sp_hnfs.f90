!!<summary>Generates the symmetry preserving HNFs for the different
!!crystal lattices.</summary>
Module sp_hnfs
  implicit none
  private
  public sc, fcc, bcc, hex, trig, st, bt

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

  !!<summary>Finds the symmetry preserving HNFs for the hexagonal
  !!lattice with determinant n. Assuming the basis of A =
  !![[1,0,0],[0.5,-0.8660254037844386,0],[0,0,2]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  SUBROUTINE hex(n,spHNFs)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)

    integer, pointer :: diagonals(:,:) => null()
    real :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, total_hnfs, status, j,k 
    integer, allocatable :: temp_HNFs(:,:,:)

    real :: beta13, beta11, gamma13, gamma11, gamma12, gamma21, gamma22
    real, allocatable :: es(:)

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0

    do i = 1,nds
       total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
    end do

    allocate(temp_HNFs(3,3,total_hnfs),STAT=status)
    if (status/=0) stop "Failed to allocate memory in hex."
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(c,a) ==0) then
          b = 0
          do while (b < c)
             beta13 = a+2*b
             beta11 = 2*b+b*b/real(a)
             if ((MOD(beta13,c)==0) .and. (MOD(beta11,c)==0)) then
                if (MOD(f,2.0)==0) then
                   allocate(es(2))
                   es = (/0.0,(f/2.0)/)
                else
                   allocate(es(1))
                   es = (/0.0/)
                end if
                do j = 1,size(es)
                   e = es(j)
                   gamma13 = (a+2*b)*e/c
                   if (MOD(gamma13,f)==0) then
                      do k=0,int(f-1)
                         d = real(k)
                         gamma11 = b*d/a -e*beta11/c
                         gamma12 = 2*d + b*d/a - e*beta11/c
                         gamma21 = c*d/c - 2*e - b*e/a
                         gamma22 = (c*d - b*e)/a
                         if ((MOD(gamma11,f)==0) .and. (MOD(gamma12,f)==0) .and. (MOD(gamma21,f)==0) .and. (MOD(gamma22,f) ==0)) then
                            nhnfs = nhnfs + 1          
                            temp_HNFs(:,:,nhnfs) = reshape((/ int(a), int(b), int(d), &
                                 0, int(c), int(e), &
                                 0, 0, int(f)/),(/3,3/))
                         end if
                      end do
                   end if
                end do
                deallocate(es)
             end if
             b = b + a
          end do
       end if
    end do

    allocate(spHNFs(3,3,nhnfs))

    spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    
  end SUBROUTINE hex

  !!<summary>Finds the symmetry preserving HNFs for the trigonal
  !!lattice with determinant n. Assuming the basis of A =
  !![[1,2,2],[2,1,2],[4,3,3]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  SUBROUTINE trig(n,spHNFs)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)

    integer, pointer :: diagonals(:,:) => null()
    real :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, total_hnfs, status, j
    integer, allocatable :: temp_HNFs(:,:,:)

    real :: beta13, beta22, beta12, gamma11, gamma12, gamma21, gamma22
    real, allocatable :: bs(:)

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0

    do i = 1,nds
       total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
    end do

    allocate(temp_HNFs(3,3,total_hnfs),STAT=status)
    if (status/=0) stop "Failed to allocate memory in trig."
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f,a)==0) then
          if (MOD(c,2.0)==0) then
             allocate(bs(2))
             bs = (/0.0, (c/2.0) /)
          else
             allocate(bs(1))
             bs = (/0.0/)
          end if

          do j = 1,size(bs)
             b = bs(j)
             beta13 = f+b*f/a
             if (MOD(beta13,c)==0) then
                e = 0.0
                do while (e <f)
                   beta22 = e + b*e/a
                   gamma21 = c + 2*e
                   gamma11 = b + 2*b*e/c
                   if ((MOD(beta22,c)==0) .and. (MOD(gamma21,f)==0) .and. (MOD(gamma11,f)==0)) then
                      d = 0.0
                      do while (d < f)
                         beta12 = -a + b + d + d*b/a
                         gamma12 = -b -d + d*d/a - e*beta12/c
                         gamma22 = -c -2*e + d*e/a - e*beta22/c
                         if ((MOD(beta12,c)==0) .and. (MOD(gamma12,f)==0) .and. (MOD(gamma22,f)==0)) then
                            nhnfs = nhnfs + 1          
                            temp_HNFs(:,:,nhnfs) = reshape((/ int(a), int(b), int(d), &
                                 0, int(c), int(e), &
                                 0, 0, int(f)/),(/3,3/))
                         end if
                         d = d + a
                      end do
                   end if
                   e = e+a
                end do
             end if
          end do
          deallocate(bs)
       end if
    end do

    allocate(spHNFs(3,3,nhnfs))

    spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    
  end SUBROUTINE trig

  !!<summary>Finds the symmetry preserving HNFs for the simple
  !!tetragonal lattice with determinant n. Assuming the basis of A =
  !![[1,0,0],[0,1,0],[0,0,2]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  SUBROUTINE st(n,spHNFs)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)

    integer, pointer :: diagonals(:,:) => null()
    real :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, total_hnfs, status, j, k, z
    integer, allocatable :: temp_HNFs(:,:,:)

    real :: beta13, gamma13, gamma12, gamma23
    real, allocatable :: bs(:), es(:)

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0

    do i = 1,nds
       total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
    end do

    allocate(temp_HNFs(3,3,total_hnfs),STAT=status)
    if (status/=0) stop "Failed to allocate memory in st."
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)
       if (MOD(c,a)==0) then
          if (MOD(c,2.0)==0) then
             allocate(bs(2))
             bs = (/0.0,(c/2.0)/)
          else
             allocate(bs(1))
             bs = (/0.0/)
          end if

          if (MOD(f,2.0)==0) then
             allocate(es(2))
             es = (/0.0,(f/2.0)/)
          else
             allocate(es(1))
             es = (/0.0/)
          end if
          do j =1,size(bs)
             b = bs(j)
             beta13 = -a + b*b/a
             if (MOD(beta13,c)==0) then
                do k = 1,size(es)
                   e = es(k)
                   gamma12 = 2*b*e/c
                   if (MOD(gamma12,f)==0) then
                      do z =0,int(f-1)
                         d = real(z)
                         gamma13 = -e*beta13/c + d*(b/a-1)
                         gamma23 = -e*(b/a+1) +d*c/a
                         if ((MOD(gamma13,f)==0) .and. (MOD(gamma23,f)==0)) then
                            nhnfs = nhnfs + 1          
                            temp_HNFs(:,:,nhnfs) = reshape((/ int(a), int(b), int(d), &
                                 0, int(c), int(e), &
                                 0, 0, int(f)/),(/3,3/))
                         end if
                      end do
                   end if
                end do
             end if
          end do
          deallocate(bs,es)
       end if
    end do

    allocate(spHNFs(3,3,nhnfs))

    spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    
  end SUBROUTINE st

  !!<summary>Finds the symmetry preserving HNFs for the body centered
  !!tetragonal lattice with determinant n. Assuming the basis of A =
  !![[-1,1,2],[1,-1,2],[1,1,-2]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  SUBROUTINE bt(n,spHNFs)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)

    integer, pointer :: diagonals(:,:) => null()
    real :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, total_hnfs, status, j, k, z, spc, size_count
    integer, allocatable :: temp_HNFs(:,:,:)

    real :: beta11, beta12, beta22, gamma11, gamma12, gamma21, gamma22
    real, allocatable :: bs(:), es(:), ds(:)

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0

    do i = 1,nds
       total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
    end do

    allocate(temp_HNFs(3,3,total_hnfs),STAT=status)
    if (status/=0) stop "Failed to allocate memory in bt."
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((a==c) .and. (a==f)) then
          allocate(bs(1),ds(1),es(1))
          bs = (/0.0/)
          ds = (/0.0/)
          es = (/0.0/)
       elseif (f==(a*c*f)) then
          allocate(bs(1))
          bs = (/0.0/)
          if (MOD(f,2.0)==0) then
             allocate(ds(2),es(2))
             ds = (/1.0, real(int(f)/2+1)/)
             es = (/1.0, real(int(f)/2+1)/)
          else
             allocate(ds(1),es(1))
             ds = (/1.0/)
             es = (/1.0/)
          end if
       else
          call smallest_prime(int(c),spc)
          size_count = int(f)/spc + 1
          allocate(ds(size_count),es(size_count))
          z = 0
          do j = 1,size_count
             ds(j) = z
             es(j) = z
             z = z + spc
          end do
          if (MOD(a,1.0)==0) then
             allocate(bs(2))
             bs = (/1.0,real(int(c)/2+1)/)
          else
             size_count = int(c)/int(a) +1
             allocate(bs(size_count))
             z = 0
             do j = 1,size_count
                bs(j) = z
                z = z + a
             end do
          end if
       end if

       if ((MOD(f,a)==0) .and. (MOD(f,c)==0)) then
          do j =1,size(bs)
             b = bs(j)
             if ((MOD((b*f),(a*c))==0) .and. (b<c)) then
                do k=1,size(es)
                   e = es(k)
                   beta22 = b*e/a
                   gamma21 = c -e*e/c
                   if ((MOD(e,c)==0) .and. (MOD(beta22,c)==0) .and. (MOD(gamma21,f)==0) .and. (MOD(e,a)==0) .and. (e<f)) then
                      do z=1,size(ds)
                         d = ds(z)
                         beta11 = -a +b +d
                         beta12 = -a +b -b*d/a
                         gamma11 = -a +b +d -e*beta11/c
                         gamma12 = -a +b +d -d*d/a -e*beta12/c
                         gamma22 = c -d*e/a + e*beta22/c
                         if ((MOD(beta11,c)==0) .and. (MOD(beta12,c)==0) .and. (MOD(gamma11,f)==0) .and. (MOD(gamma12,f)==0) .and. (MOD(gamma22,f)==0) .and. (MOD(d,a)==0) .and. (d<f)) then
                            nhnfs = nhnfs + 1          
                            temp_HNFs(:,:,nhnfs) = reshape((/ int(a), int(b), int(d), &
                                 0, int(c), int(e), &
                                 0, 0, int(f)/),(/3,3/))
                            
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
       deallocate(bs,es,ds)
    end do

    allocate(spHNFs(3,3,nhnfs))

    spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    
  end SUBROUTINE bt





  
  !!<summary>Finds the smallest prime factor of the given positive integer.</summary>
  !!<parameter name="a" regular="true">A positive integer number.</parameter>
  !!<parameter name="sp" regular="true">The smallest prime factor of a.</parameter>
  SUBROUTINE smallest_prime(a,sp)
    integer, intent(in) :: a
    integer, intent(out) :: sp

    integer ::i

    if (a<0) stop "smallest_prime is only designed for positive integers."
    if (a <= 2) then
       sp = a
    else
       do i = 2, a
          if (MOD(a,i)==0) then
             sp = i
             exit
          end if
       end do
    end if

  end SUBROUTINE smallest_prime

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
    if(status/=0) stop "Allocation failed in get_HNF_diagonals"
    diagonals = tempDiag(:,1:id)
  END SUBROUTINE get_HNF_diagonals
  
end Module sp_hnfs
