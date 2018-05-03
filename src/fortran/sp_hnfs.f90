!!<summary>Generates the symmetry preserving HNFs for the different
!!crystal lattices.</summary>
Module sp_hnfs
  use grid_utils
  
  implicit none
  private
  public sc_3, fcc_1, bcc_5, hex_12, hex_22, rhom_9, rhom_4_2, rhom_24, st_11, st_21, &
       bct_15, bct_6, bct_7, bct_18, so_32, baseco_23, baseco_36, baseco_40, baseco_38_13, &
       bco_19, bco_8, bco_42, fco_26, fco_16, sm_33, sm_34, sm_35, basecm_43, basecm_14, &
       basecm_28, basecm_41, basecm_27, basecm_37_39, basecm_29_30, basecm_10_17, &
       basecm_20_25, tric_31_44

  integer, parameter:: dp=selected_real_kind(15,307)
  integer, parameter:: sp=selected_real_kind(6,37)
  integer, parameter:: si=selected_int_kind(1) ! very short integer -10..10 range
  integer, parameter:: li=selected_int_kind(18) ! Big integer -10^18..10^18 range
  
CONTAINS

  !!<summary>Finds the symmetry preserving HNFs for the simple cubic
  !!lattice with determinant n. Assuming the basis of A =
  !![[1,0,0],[0,1,0],[0,0,1]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grids
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  SUBROUTINE sc_3(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, n_irr, nhnfs, eps_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable, intent(out) :: grids(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i
    integer :: temp_HNFs(3,3,1)
    real(dp) :: eps

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if
    
    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)

    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((a==c) .and. (c==f)) then
          b = 0
          d = 0
          e = 0
          temp_HNFs(:,:,1) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((a==c) .and. (real(f,dp)/real(a,dp)==2)) then
          b = 0
          d = a
          e = a
          temp_HNFs(:,:,1) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((c==f) .and. (real(f,dp)/real(a,dp)==2)) then
          b = a
          d = a
          e = 0
          temp_HNFs(:,:,1) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       end if
    end do
    allocate(spHNFs(3,3,1))
    nhnfs = 1
    spHNFs(:,:,1) = temp_HNFs(:,:,1)
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp
    call grid_metrics(U, B_vecs, at, spHNFs(:,:,1), No, Nu, Co, Cu, O, grid, rmin, n_irr, eps_=eps)
    
  end SUBROUTINE sc_3

  !!<summary>Finds the symmetry preserving HNFs for the face centered
  !!cubic lattice with determinant n. Assuming the basis of A =
  !![[0,1,1],[1,0,1],[1,1,0]].</summary>  
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grids
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  SUBROUTINE fcc_1(n, No, Nu, Co, Cu, O, U, B_vecs, spHNFs, grid, rmin, n_irr, eps_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable, intent(out) :: grids(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i
    integer :: temp_HNFs(3,3,1)
    real(dp) :: eps

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if
    
    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)

    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((a==c) .and. (a==f)) then
          b = 0
          d = 0
          e = 0
          temp_HNFs(:,:,1) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((c==f) .and. ((real(f,dp)/real(a,dp) == 2) .or. (real(f,dp)/real(a,dp) == 4))) then
          b = a
          d = a
          e = 0
          temp_HNFs(:,:,1) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       end if
    end do
    allocate(spHNFs(3,3,1))
    nhnfs = 1    
    spHNFs(:,:,1) = temp_HNFs(:,:,1)
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp
    allocate(grids(3,3,1),STAT=status)
    call grid_metrics(U, B_vecs, at, spHNFs(:,:,1), No, Nu, Co, Cu, O, grid, rmin, n_irr, eps_=eps)
    
  end SUBROUTINE fcc_1

  !!<summary>Finds the symmetry preserving HNFs for the body centered
  !!cubic lattice with determinant n. Assuming the basis of A =
  !![[-1,1,1],[1,-1,1],[1,1,-1]].</summary>  
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grids
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  SUBROUTINE bcc_5(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, n_irr, nhnfs, eps_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i
    integer :: temp_HNFs(3,3,1)
    real(dp) :: eps

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if
    
    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)

    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((a==c) .and. (a==f)) then
          b = 0
          d = 0
          e = 0
          temp_HNFs(:,:,1) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((a==c) .and. (real(f,dp)/real(a,dp)==2)) then
          b = 0
          d = a
          e = a
          temp_HNFs(:,:,1) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((a==c) .and. (real(f,dp)/real(a,dp) ==4)) then
          b = 0
          d = 3*a
          e = 3*a
          nhnfs = nhnfs + 1          
          temp_HNFs(:,:,1) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       end if
    end do
    allocate(spHNFs(3,3,nhnfs))
    nhnfs = 1
    spHNFs(:,:,1) = temp_HNFs(:,:,1)
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp
    call grid_metrics(U, B_vecs, at, spHNFs(:,:,1), No, Nu, Co, Cu, O, grid, rmin, n_irr, eps_=eps)
    
  end SUBROUTINE bcc_5

  !!<summary>Finds the symmetry preserving HNFs for the hexagonal
  !!lattice with determinant n. Assuming the basis of A =
  !![[1,0,0],[0.5,-0.8660254037844386,0],[0,0,2]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE hex_12(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, status, j,k, nes
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta13, beta11, gamma13, gamma11, gamma12, gamma21, gamma22
    integer :: es(2)
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(c, a) ==0) then
          b = 0
          do while (b < c)
             if (MOD(b*b, a)==0) then
                beta13 = a+2*b
                beta11 = 2*b+b*b/a
                if ((MOD(beta13, c)==0) .and. (MOD(beta11, c)==0)) then
                   if (MOD(f, 2)==0) then
                      es = (/0, (f/2)/)
                      nes = 2
                   else
                      es = (/0, -1/)
                      nes = 1
                   end if
                   do j = 1, nes
                      e = es(j)
                      if (MOD((a+2*b)*e, c)==0) then
                         gamma13 = (a+2*b)*e/c
                         if (MOD(gamma13, f)==0) then
                            do k=0, (f-1)
                               d = k
                               if ((MOD(b*d, a)==0) .and. (MOD(e*beta11, c)==0) .and. &
                                    (MOD(b*e, a)==0) .and. (MOD(c*d, a)==0) .and. &
                                    (MOD(c*d-b*e, a)==0)) then 
                                  gamma11 = b*d/a -e*beta11/c
                                  gamma12 = 2*d + b*d/a - e*beta11/c
                                  gamma21 = c*d/a - 2*e - b*e/a
                                  gamma22 = (c*d - b*e)/a
                                  if ((MOD(gamma11,f)==0) .and. (MOD(gamma12,f)==0) .and. &
                                       (MOD(gamma21,f)==0) .and. (MOD(gamma22,f) ==0)) then
                                     nhnfs = nhnfs + 1
                                     if (all_hnfs) then 
                                        temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                     else
                                        
                                        temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                        call compare_grids(lat_vecs, B_vecs, at, &
                                             temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                             grid, rmin, n_irr, eps_)
                                     end if
                                  end if
                               end if
                            end do
                         end if
                      end if
                   end do
                end if
                b = b + a
             end if
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
  end SUBROUTINE hex_12

  !!<summary>Finds the symmetry preserving HNFs for the hexagonal
  !!lattice with determinant n. Assuming the basis of A =
  !![[0,0,-0.5],[1,0,0],[-0.5,0.8660254037844386,0]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE hex_22(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma21, gamma22, gamma11, gamma12
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_22."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, c) ==0) then
          do e=0, (f-1), c
             if (MOD(e*e, c)==0) then
                gamma21 = -c+(e*2)
                gamma22 = -c+e-(e*e)/c
                if ((MOD(gamma21, f)==0) .and. (MOD(gamma22, f)==0)) then
                   do b=0, (c-1)
                      do d=0, (f-1), c
                         if (MOD(d*e, c)==0) then
                            gamma11 = -b+d*2
                            gamma12 = -b+d-(d*e)/c
                            if ((MOD(gamma11, f)==0) .and. (MOD(gamma12, f)==0)) then
                               nhnfs = nhnfs + 1          
                               if (all_hnfs) then 
                                  temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                               else
                                  
                                  temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                                  call compare_grids(lat_vecs, B_vecs, at, &
                                       temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                       grid, rmin, n_irr, eps_)
                               end if
                            end if
                         end if
                      end do
                   end do
                end if
             end if
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
  end SUBROUTINE hex_22

  !!<summary>Finds the symmetry preserving HNFs for the rhombohedral
  !!lattice with determinant n. Assuming the basis of A =
  !![[1,2,2],[2,1,2],[4,3,3]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE rhom_9(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta13, beta22, beta12, gamma11, gamma12, gamma21, gamma22
    integer :: bs(2), nbs
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in rhom_9."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f,a)==0) then
          if (MOD(c,2)==0) then
             bs = (/0, (c/2) /)
             nbs = 2
          else
             bs = (/0, -1/)
             nbs = 1
          end if

          do j = 1, nbs
             b = bs(j)
             if (MOD(b*f, a)==0) then
                beta13 = f+b*f/a
                if (MOD(beta13,c)==0) then
                   e = 0
                   do while (e <f)
                      if ((MOD(b*e, a)==0) .and. (MOD(2*b*e, c)==0)) then
                         beta22 = e + b*e/a
                         gamma21 = c + 2*e
                         gamma11 = b + 2*b*e/c
                         if ((MOD(beta22, c)==0) .and. (MOD(gamma21, f)==0) .and. &
                              (MOD(gamma11, f)==0)) then
                            d = 0
                            do while (d < f)
                               if ((MOD(d*b, a)==0) .and. (MOD(d*d, a)==0) .and. &
                                    (MOD(e*beta12, c)==0) .and. (MOD(d*e, a)==0) &
                                    .and. (MOD(d*e, a)==0) .and. (MOD(e*beta22, c)==0)) then
                                  beta12 = -a + b + d + d*b/a
                                  gamma12 = -b -d + (d*d/a) - e*beta12/c
                                  gamma22 = -c -2*e + (d*e/a) - e*beta22/c
                                  if ((MOD(beta12, c)==0) .and. (MOD(gamma12, f)==0) .and. &
                                       (MOD(gamma22, f)==0)) then
                                     nhnfs = nhnfs + 1          
                                     if (all_hnfs) then 
                                        temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                     else
                                        
                                        temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                        call compare_grids(lat_vecs, B_vecs, at, &
                                             temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                             grid, rmin, n_irr, eps_)
                                     end if
                                  end if
                               end if
                               d = d + a
                            end do
                         end if
                      end if
                      e = e+a
                   end do
                end if
             end if
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE rhom_9

  !!<summary>Finds the symmetry preserving HNFs for the rhombohedral
  !!lattice with determinant n. Assuming A = [[-1, 0,-1],[0, -1.32288,
  !!-0.5],[-1.11652, -0.610985, 0.616515]] for number 2 and A= [[-1,
  !!0,-1],[0, -1.32288, 0.5],[-0.548584, 0.774292, 1.04858]] for
  !!number 4.</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE rhom_4_2(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta32, beta22, beta11, beta12, gamma11, gamma12, gamma22
    integer :: nbs, nds
    integer, allocatable :: bs(:), ds(:)
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in rhom_4_2."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((MOD(f, a)==0) .and. (MOD(f, c)==0))then
          if (c<f) then
             e = f-c
          else
             e = 0
          end if
          allocate(bs(c))
          do j=0, (c-1)
             bs(j+1) = j
             nbs = c
          end do
          if (.not. e==0) then
             if (a<c) then
                bs = -1
                nbs = 2
                bs(1) = 0
                bs(2) = a
             else
                bs = -1
                nbs = 1
                bs(1) = 0
             end if
          else
             bs = -1
             nbs = 1
             if (a<c) then
                bs(1) = a
             else
                bs(1) = 0
             end if
          end if
          if (c==1) then
             allocate(ds(1))
             nds = 1
             ds = e
          else if (e==0) then
             allocate(ds(1))
             nds = 1
             if (a==c) then
                ds = 0
             else 
                ds = a
             end if
          else 
             allocate(ds(f))
             nds = f
             do j=1, f
                ds(j) = j-1
             end do
          end if

          if (MOD(e*e, c)==0) then
             gamma12 = -c+e*e/c
             if ((MOD(e, c)==0) .and. (MOD(gamma12, f)==0) .and. (MOD(e, a)==0)) then
                do j=1, nbs
                   b = bs(j)
                   if ((MOD(b*f, a)==0) .and. (MOD(b*e, a)==0)) then
                      beta32 = b*f/a
                      beta22 = b*e/a
                      if ((MOD(beta32, c)==0) .and. (MOD(beta22, c)==0)) then
                         do k=1, nds
                            d = ds(k)
                            if ((MOD(b*d, a)==0) .and. (MOD(beta11*e, c)==0) .and. &
                                 (MOD(d*d, a)==0) .and. (MOD(beta12*e, c)==0) .and. &
                                 (MOD(d*e, a)==0) .and. (MOD(e*e, a*c)==0)) then
                               beta11 = b-d
                               beta12 = -a+b*d/a
                               gamma11 = -b+d-beta11*e/c
                               gamma12 = -b+(d*d/a)-beta12*e/c
                               gamma22 = -c+(d*e/a)-b*e*e/(a*c)
                               if ((MOD(beta11, c)==0) .and. (MOD(beta12, c)==0) .and. &
                                    (MOD(gamma11, f)==0) .and. (MOD(gamma12, f)==0) .and. &
                                    (MOD(gamma22, f)==0)) then 
                                  nhnfs = nhnfs + 1          
                                  if (all_hnfs) then 
                                     temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                          0, c, e, &
                                          0, 0, f/),(/3,3/))
                                  else
                                     
                                     temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                          0, c, e, &
                                          0, 0, f/),(/3,3/))
                                     call compare_grids(lat_vecs, B_vecs, at, &
                                          temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                          grid, rmin, n_irr, eps_)
                                  end if
                               end if
                            end if
                         end do
                      end if
                   end if
                end do
             end if
          end if
          deallocate(bs,ds)
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE rhom_4_2


  !!<summary>Finds the symmetry preserving HNFs for the rhombohedral
  !!lattice with determinant n. Assuming A =
  !![[-1,0,-1],[1.51184,0,-0.845178],[-0.255922,-1.44338,0.92259]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE rhom_24(n,spHNFs)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma11, gamma21, gamma22
    integer, allocatable :: ds(:), temp(:)
    integer :: count
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in rhom_24."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       b = 0
       
       if ((MOD(f, a)==0) .and. (MOD(f, c)==0)) then
          if (c==f) then
             allocate(ds(1))
             e = 0
             ds = 0
          else
             e = f-c
             if ((c==1) .and. (f>3)) then
                allocate(ds(1))
                ds = 3
             else
                allocate(temp(f))
                count = 0
                do j=1, f
                   if ((MOD(j-1, c)==0) .and. ((j-1)<(e+1))) then
                      count = count + 1
                      temp(count) = j-1
                   end if
                end do
                allocate(ds(count))
                ds(1:count) = temp(1:count)
                deallocate(temp)                
             end if
          end if

          if ((MOD(e, c)==0) .and. (MOD(e, a)==0)) then
             do k=1,size(ds)
                d = ds(k)
                if ((MOD(d*e, c)==0) .and. (MOD(d*d, a)==0) .and. (MOD(e*e, c) == 0) &
                     .and (MOD(d*e, a)==0)) then
                   gamma11 = 2*d-(d*d/a)-d*e/c
                   gamma21 = 2*e-(d*e/a)-e*e/c
                   gamma22 = -c+e-(d*e/a)-e*e/c
                   if ((MOD(gamma11, f)==0) .and. (MOD(gamma21, f)==0) .and. &
                        (MOD(gamma22, f)==0)) then
                      nhnfs = nhnfs + 1          
                      if (all_hnfs) then 
                         temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                      else
                                     
                         temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                         call compare_grids(lat_vecs, B_vecs, at, &
                              temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                              grid, rmin, n_irr, eps_)
                      end if
                   end if
                end if
             end do
          end if
          deallocate(ds)
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE rhom_24
  
  !!<summary>Finds the symmetry preserving HNFs for the simple
  !!tetragonal lattice with determinant n. Assuming the basis of A =
  !![[1,0,0],[0,1,0],[0,0,2]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE st_11(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta13, gamma13, gamma12, gamma23
    integer :: bs(2), es(2), nbs, nes
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in st_11."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)
       
       if (MOD(c, a)==0) then
          if (MOD(c, 2)==0) then
             bs = (/0, (c/2)/)
             nbs = 2
          else
             bs = (/0, -1/)
             nbs = 1
          end if

          if (MOD(f, 2)==0) then
             nes = 2
             es = (/0, (f/2)/)
          else
             nes = 1
             es = (/0, -1/)
          end if
          do j =1, nbs
             b = bs(j)
             if (MOD(b*b, a)==0) then
                beta13 = -a + b*b/a
                if (MOD(beta13, c)==0) then
                   do k = 1, nes
                      e = es(k)
                      if (MOD(2*b*e, c)==0) then
                         gamma12 = 2*b*e/c
                         if (MOD(gamma12, f)==0) then
                            do d = 0, (f-1)
                               if ((MOD(b, a)==0) .and. (MOD(d*c, a)==0)) then
                                  gamma13 = -e*beta13/c + d*(b/a-1)
                                  gamma23 = -e*(b/a+1) +d*c/a
                                  if ((MOD(gamma13,f)==0) .and. (MOD(gamma23,f)==0)) then
                                     nhnfs = nhnfs + 1          
                                     if (all_hnfs) then 
                                        temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                     else
                                     
                                        temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                        call compare_grids(lat_vecs, B_vecs, at, &
                                             temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                             grid, rmin, n_irr, eps_)
                                     end if
                                  end if
                               end if
                            end do
                         end if
                      end if
                   end do
                end if
             end if
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE st_11
  
  !!<summary>Finds the symmetry preserving HNFs for the simple
  !!tetragonal lattice with determinant n. Assuming the basis of A =
  !![[0,0,0.5],[1,0,0],[0,1,0]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE st_21(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta13, gamma13, gamma12, gamma23
    integer :: bs(2), es(2), nbs, nes
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in st_21."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, c)==0) then
          if (MOD(c, 2)==0) then
             bs = (/0, (c/2)/)
             nbs = 2
          else
             bs = (/0, -1/)
             nbs = 2
          end if
          do j=1, nbs
             b = bs(j)
             if (MOD(f, 2)==0) then
                es = (/0, (f/2)/)
                nes = 2
             else
                es = (/0, -1/)
                nes = 1
             end if
             do k=1, nes
                d = es(k)
                beta13 = b-d
                if (MOD(beta13, c)==0) then
                   do z=1, nes
                      e = es(z)
                      if ((MOD(2*d*e, c)==0) .and. (MOD(e*e, c)==0)) then
                         gamma12 = 2*d-2*d*e/c
                         gamma13 = -b+d-e*beta13/c
                         gamma23 = -c+e*e/c
                         if ((MOD(gamma12, f)==0) .and. (MOD(gamma13, f)==0) .and. &
                              (MOD(gamma23, f)==0)) then
                            nhnfs = nhnfs + 1          
                            if (all_hnfs) then 
                               temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                            else
                               
                               temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                               call compare_grids(lat_vecs, B_vecs, at, &
                                    temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                    grid, rmin, n_irr, eps_)
                            end if
                         end if
                      end if
                   end do
                end if
             end do
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE st_21
  
  !!<summary>Finds the symmetry preserving HNFs for the body centered
  !!tetragonal lattice with determinant n. Assuming the basis of A =
  !![[-1,1,2],[1,-1,2],[1,1,-2]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE bct_15(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z, spc, size_count
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta11, beta12, beta22, gamma11, gamma12, gamma21, gamma22
    integer, allocatable :: bs(:), es(:), ds(:)
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in bct_15."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((a==c) .and. (a==f)) then
          allocate(bs(1), ds(1), es(1))
          bs = (/0/)
          ds = (/0/)
          es = (/0/)
       elseif (f==(a*c*f)) then
          allocate(bs(1))
          bs = (/0/)
          if (MOD(f, 2)==0) then
             allocate(ds(2),es(2))
             ds = (/1, (f/2)+1/)
             es = (/1, (f/2)+1/)
          else
             allocate(ds(1), es(1))
             ds = (/1/)
             es = (/1/)
          end if
       else
          call smallest_prime(c, spc)
          size_count = f/spc + 1
          allocate(ds(size_count), es(size_count))
          z = 0
          do j = 1, size_count
             ds(j) = z
             es(j) = z
             z = z + spc
          end do
          if (a==1) then
             allocate(bs(2))
             bs = (/1, (c/2)+1/)
          else
             size_count = c/a +1
             allocate(bs(size_count))
             z = 0
             do j = 1, size_count
                bs(j) = z
                z = z + a
             end do
          end if
       end if

       if ((MOD(f, a)==0) .and. (MOD(f, c)==0)) then
          do j =1, size(bs)
             b = bs(j)
             if ((MOD((b*f), (a*c))==0) .and. (b<c)) then
                do k=1, size(es)
                   e = es(k)
                   if ((MOD(b*e, a)==0) .and. (MOD(e*e, c)==0)) then
                      beta22 = b*e/a
                      gamma21 = c -e*e/c
                      if ((MOD(e, c)==0) .and. (MOD(beta22, c)==0) .and. &
                           (MOD(gamma21, f)==0) .and. (MOD(e, a)==0) .and. (e<f)) then
                         do z=1, size(ds)
                            d = ds(z)
                            if ((MOD(b*d, a)==0) .and. (MOD(e*beta11, c)==0) .and. &
                                 (MOD(d*d, a)==0) .and. (MOD(e*beta12, c)==0)) then
                               beta11 = -a +b +d
                               beta12 = -a +b -b*d/a
                               gamma11 = -a +b +d -e*beta11/c
                               gamma12 = -a +b +d -d*d/a -e*beta12/c
                               gamma22 = c -d*e/a + e*beta22/c
                               if ((MOD(beta11, c)==0) .and. (MOD(beta12, c)==0) .and. &
                                    (MOD(gamma11, f)==0) .and. (MOD(gamma12, f)==0) .and. &
                                    (MOD(gamma22, f)==0) .and. (MOD(d, a)==0) .and. (d<f)) then
                                  nhnfs = nhnfs + 1          
                                  if (all_hnfs) then 
                                     temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                          0, c, e, &
                                          0, 0, f/),(/3,3/))
                                  else
                                     
                                     temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                          0, c, e, &
                                          0, 0, f/),(/3,3/))
                                     call compare_grids(lat_vecs, B_vecs, at, &
                                          temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                          grid, rmin, n_irr, eps_)
                                  end if
                               end if
                            end if
                         end do
                      end if
                   end if
                end do
             end if
          end do
       end if
       deallocate(bs,es,ds)
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE bct_15
  
  !!<summary>Finds the symmetry preserving HNFs for the body centered
  !!tetragonal lattice with determinant n. Assuming the basis of A =
  !![[-1.95095 , 1.41625 , -0.433603], [ 1.  , -1.  , -2.  ],
  !![ 1.95095 , 1.19163 , 0.879663]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE bct_7(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta11, beta33, beta23, alpha23, alpha13, beta13
    integer :: gamma13, gamma23, gamma11
    integer, allocatable :: bs(:), es(:), temp(:), ds(:)
    integer :: count
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in bct_7."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((MOD(f, a)==0) .and. (MOD(f, c)==0) .and. (MOD(a, c)==0)) then
          if ((a==1) .and. (c==1)) then
             allocate(es(1), bs(1))
             es = (/f-1/)
             bs = (/0/)
          else if ((a==c) .and. (c==f)) then
             allocate(es(1), bs(1))
             es = (/0.0_dp/)
             bs = (/0.0_dp/)
          else
             allocate(bs(c), temp(f))
             count = 0
             do j=0, f-1
                if (j==0) then
                   count = count + 1
                   temp(count) = j
                   count = count + 1
                   temp(count) = c
                else if ((c+a*j)<f) then
                   count = count + 1
                   temp(count) = c+a*j
                end if
             end do
             allocate(es(count))
             es(1:count) = temp(1:count)
             deallocate(temp)
             do j=0, (c-1)
                bs(j+1) = j
             end do
          end if

          do j=1, size(bs)
             b = bs(j)
             if (MOD(b*f, a)==0) then
                beta11 = -a+2*b
                beta33 = f-b*f/a
                if ((MOD(beta11, c)==0) .and. (MOD(beta33, c)==0)) then
                   do k=1, size(es)
                      e = es(k)
                      if (MOD(alpha23*b, c)==0) then
                         alpha23 = -c+e
                         beta23 = e-b*alpha23/a
                         if ((MOD(alpha23, a)==0) .and. (MOD(beta23, c)==0)) then
                            if (b==0) then
                               if ((MOD(f, 2)==0) .and. (((f/2)+a)<f)) then
                                  allocate(ds(3))
                                  ds = (/0, a, ((f/2)+a)/)
                               else if (a<f) then
                                  allocate(ds(2))
                                  ds = (/0, a/)
                               else
                                  allocate(ds(int(f)))
                                  do z=0, (f-1)
                                     ds(z+1) = z
                                  end do
                               end if
                            else
                               if ((e==(f-c)) .and. (MOD(f, 2)==0)) then
                                  allocate(ds(1))
                                  ds = (/(f/2)+b/)
                               else 
                                  allocate(ds(f))
                                  do z=0, (f-1)
                                     ds(z+1) = z
                                  end do
                               end if
                            end if
                            do z=1, size(ds)
                               d = ds(z)
                               alpha13 = -b+d
                               if (MOD(b*alpha13, a)==0) then
                                  beta13 = -a+d-b*alpha13/a
                                  if (MOD(beta13*e, c)==0) then
                                     gamma11 = -a+2*d-e*beta11/c
                                     gamma13 = d-(d*alpha13/a)-e*beta13/c
                                     gamma23 = e-(d*alpha23/a)-e*beta23/c
                                     if ((MOD(alpha13, a)==0) .and. (MOD(beta13, c)==0) &
                                          .and. (MOD(gamma11, f)==0) .and. &
                                          (MOD(gamma13, f)==0) .and. (MOD(gamma23, f)==0)) then
                                        nhnfs = nhnfs + 1          
                                        if (all_hnfs) then 
                                           temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                                0, c, e, &
                                                0, 0, f/),(/3,3/))
                                        else
                                           
                                           temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                                0, c, e, &
                                                0, 0, f/),(/3,3/))
                                           call compare_grids(lat_vecs, B_vecs, at, &
                                                temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                                grid, rmin, n_irr, eps_)
                                        end if
                                     end if
                                  end if
                               end if
                            end do
                            deallocate(ds)
                         end if
                      end if
                   end do
                end if
             end if
          end do
          deallocate(es,bs)
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE bct_7  
  
  !!<summary>Finds the symmetry preserving HNFs for the body centered
  !!tetragonal lattice with determinant n. Assuming the basis of A =
  !![[-1.  , 1.  , 2.  ], [ 1.  , 1.60788 , -1.55394 ], [ 1.80278 ,
  !!-1.47253 , 0.762655]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE bct_6(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta11, beta12, beta22, gamma11, gamma12, gamma21, gamma22
    integer, allocatable :: bs(:), es(:), ds(:)
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in bct_6."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((MOD(f, c)==0) .and. (MOD(f, a)==0)) then
          if ((a==c) .and. (f==a)) then
             allocate(es(1), bs(1), ds(1))
             es = (/0/)
             ds = (/0/)
             bs = (/0/)
          else
             allocate(bs(c))
             do j=0, (c-1)
                bs(j+1) = j
             end do
             if (((f/2)+c)<f) then
                allocate(es(3))
                es = (/0, c, ((f/2)+c)/)
             else if (c<f) then 
                allocate(es(2))
                es = (/0, c/)
             else
                allocate(es(1))
                es = (/0/)
             end if
             if (MOD(c, 2)==0) then
                allocate(ds((f/(c/2))))
                do j=1, (f/(c/2))
                   ds(j) = (j-1)*(c/2)
                end do
             else
                allocate(ds(f/c))
                do j=1, (f/c)
                   ds(j) = (j-1)*c
                end do
             end if
          end if

          do j=1, size(bs)
             b = bs(j)
             if (MOD(b*f, a*c)==0) then
                do k=1, size(es)
                   e = es(k)
                   if ((MOD(e*e, c)==0) .and. (MOD(b*e, a)==0))then
                      gamma21 = c-e*e/c
                      beta22 = b*e/a
                      if ((MOD(e, a)==0) .and. (MOD(gamma21, f)==0) .and. &
                           (MOD(beta22, c)==0)) then
                         do z=1, size(ds)
                            d = ds(z)
                            beta11 = -a+b+d
                            if ((MOD(b*d, a)==0) .and.(MOD(d*d, a)==0) .and. &
                                 (MOD(d*e, a)==0) .and. (MOD(beta11*e, c)==0)) then
                               beta12 = -a+b-b*d/a
                               if (MOD(beta12*e, c==0)) then
                                  gamma11 = beta11-beta11*e/c
                                  gamma12 = beta11 -(d*d/a)-beta12*e/c
                                  gamma22 = c-(d*e/a)+e*beta22/c
                                  if ((MOD(beta11, c)==0) .and. (MOD(beta12, c)==0) .and. &
                                       (MOD(gamma11, f)==0) .and. (MOD(gamma12, f)==0) .and. &
                                       (MOD(gamma22, f)==0)) then
                                     nhnfs = nhnfs + 1          
                                     if (all_hnfs) then 
                                        temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                     else
                                        
                                        temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                        call compare_grids(lat_vecs, B_vecs, at, &
                                             temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                             grid, rmin, n_irr, eps_)
                                     end if
                                  end if
                               end if
                            end if
                         end do
                      end if
                   end if
                end do
             end if
          end do
          deallocate(es,bs,ds)
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE bct_6
  
  !!<summary>Finds the symmetry preserving HNFs for the body centered
  !!tetragonal lattice with determinant n. Assuming the basis of A =
  !![[ 0, 0, 2], [ 1, -2, 1], [-2, -1, 1]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE bct_18(n,spHNFs)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z, sp
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta31, beta33, beta23, beta21, beta22, beta11, beta12, beta13
    integer :: gamma11, gamma12, gamma21, gamma23, gamma13, alpha11, alpha21, gamma22
    integer, allocatable :: bs(:), es(:), ds(:)
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in bct_18."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if ((MOD(c, a)==0) .and. (MOD(f, a)==0) .and. (MOD(f, c)==0)) then
          if (c==f) then
             allocate(es(1))
             es = (/0/)
             if (a==c) then
                allocate(bs(1), ds(1))
                bs = (/0/)
                ds = (/0/)
             else
                call smallest_prime(c,sp)
                allocate(bs(c/sp), ds(f/sp))
                do j=1, size(bs)
                   bs(j) = (j-1)*sp
                end do
                do j=1, size(ds)
                   ds(j) = (j-1)*sp,dp
                end do
             end if
          else if ((a*c*f)==f) then
             allocate(bs(1))
             bs = (/0/)
             if (MOD(f, 2)==0) then
                allocate(es(2))
                es = (/(f/2)-1, f-1/)
             else
                allocate(es(1))
                es = (/f-1/)
             end if
             allocate(ds(1))
             if (f>=2) then
                ds = (/f-2/)
             else
                ds = (/0/)
             end if
          else
             call smallest_prime(c, sp)
             allocate(bs(c/sp))
             do j=1, size(bs)
                bs(j) = (j-1)*sp
             end do
             if (MOD(f, 2)==0) then
                allocate(es(2))
                es = (/(f/2)-c, f-c/)
             else
                allocate(es(1))
                es = (/f-c/)
             end if
             if (MOD(c, 2)==0) then
                allocate(ds(f/(c/2)))
                do j=1, size(ds)
                   ds(j) = (j-1)*(c/2)
                end do
             else
                allocate(ds(f/c))
                do j=1, size(ds)
                   ds(j) = (j-1)*c
                end do
             end if
          end if

          do j=1, size(bs)
             b = bs(j)
             if (MOD(b*f, a)==0) then
                beta31 = f+b*f/a
                beta33 = b*f/a
                if ((MOD(beta31, c)==0) .and. (MOD(beta33, c)==0)) then
                   do k=1, size(es)
                      e = es(k)
                      alpha21 = -c-e
                      if ((MOD(b*alpha21, a)==0) .and. (MOD(b*c, a)==0)) then
                         beta23 = -b*alpha21/a
                         beta21 = beta23+e
                         beta22 = (b*c/a)-e
                         if ((MOD(alpha21, a)==0) .and. (MOD(beta23, c)==0) .and. &
                              (MOD(beta21, c)==0) .and. (MOD(beta22, c)==0) .and. &
                              (MOD(beta23, c)==0)) then
                            do z=1,size(ds)
                               d = ds(z)
                               alpha11 = -b-d
                               if ((MOD(b*alpha11, a)==0) .and. (MOD(b*b, a)==0) .and. &
                                    (MOD(alpha11*d, a)==0) .and. (MOD(b*d, a)==0) .and. &
                                    (MOD(c*d, a)==0)) then
                                  beta11 = b-(b*alpha11/a)+d
                                  beta12 = b+(b*b/a)-d
                                  beta13 = 2*b-b*alpha11/a
                                  if ((MOD(beta11*e, c)==0) .and. (MOD(beta13*e, c)==0) .and. &
                                       (MOD(beta12*e, c)==0)) then
                                     gamma11 = b+d-(alpha11*d/a)-beta11*e/c
                                     gamma12 = b+d+(b*d/a)-beta12*e/c
                                     gamma13 = 2*d-(alpha11*d/a)-beta13*e/c
                                     gamma21 = c-(d*alpha21/a)-e*beta21/c
                                     gamma22 = c+(c*d/a)-beta22*e/c
                                     gamma23 = -(d*alpha21/a)-beta23*e/c
                                     if ((MOD(alpha11,a)==0) .and. (MOD(beta11,c)==0) .and. &
                                          (MOD(beta12,c)==0) .and. (MOD(beta13,c)==0) .and. &
                                          (MOD(gamma11,f)==0) .and. (MOD(gamma12,f)==0) .and. &
                                          (MOD(gamma13,f)==0)  .and. (MOD(gamma21,f)==0) .and. &
                                          (MOD(gamma22,f)==0) .and. (MOD(gamma23,f)==0)) then
                                        nhnfs = nhnfs + 1          
                                        if (all_hnfs) then 
                                           temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                                0, c, e, &
                                                0, 0, f/),(/3,3/))
                                        else
                                           
                                           temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                                0, c, e, &
                                                0, 0, f/),(/3,3/))
                                           call compare_grids(lat_vecs, B_vecs, at, &
                                                temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                                grid, rmin, n_irr, eps_)
                                        end if
                                     end if
                                  end if
                               end if
                            end do
                         end if
                      end if
                   end do
                end if
             end if
          end do
          deallocate(es,ds,bs)
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE bct_18
  
  !!<summary>Finds the symmetry preserving HNFs for the simple
  !!orthorhombic lattice with determinant n. Assuming the basis of A =
  !![[1,0,0],[0,2,0],[0,0,3]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE so_32(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: bs(2), es(2), ds(2), nbs, ne_ds
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(c, 2)==0) then
          nbs = 2
          bs = (/0, c/2/)
       else
          nbs = 1
          bs = (/0, -1/)
       end if
       if (MOD(f, 2)==0) then
          ne_ds = 2
          es = (/0, f/2/)
          ds = (/0, f/2/)
       else
          ne_ds = 1
          es = (/0, -1/)
          ds = (/0, -1/)          
       end if

       do j = 1, nbs
          b = bs(j)
          do k=1, ne_ds
             e = es(k)
             if (MOD((2*b*e),(f*c))==0) then
                do z = 1,size(ds)
                   d = ds(z)
                   nhnfs = nhnfs + 1          
                   if (all_hnfs) then 
                      temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                           0, c, e, &
                           0, 0, f/),(/3,3/))
                   else
                      
                      temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                           0, c, e, &
                           0, 0, f/),(/3,3/))
                      call compare_grids(lat_vecs, B_vecs, at, &
                           temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                           grid, rmin, n_irr, eps_)
                   end if
                end do
             end if
          end do
       end do
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE so_32

  !!<summary>Finds the symmetry preserving HNFs for the face centered
  !!orthorhombic lattice with determinant n. Assuming the basis of A =
  !![[0,1,1.5],[0.5,0,1.5],[0,0,3]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE fco_26(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma12, gamma13, gamma23
    integer :: bs(2), nbs
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(c, 2)==0) then
          nbs = 2
          bs = (/0, (c/2)/)
       else
          nbs = 1
          bs = (/0, -1/)
       end if

       do j=1, nbs
          b = bs(j)
          do e = 0, (f-1)
             if (MOD(2*b*e, c) == 0) then
                gamma23 = c +2*e
                gamma12 = b + 2*b*e/c
                if ((MOD(gamma23, f)==0) .and. (MOD(gamma12, f)==0)) then
                   do d = 0, (f-1)
                      gamma13 = a +b +2*d
                      if (MOD(gamma13, f)==0) then
                         nhnfs = nhnfs + 1          
                         if (all_hnfs) then 
                            temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                         else
                                     
                            temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                            call compare_grids(lat_vecs, B_vecs, at, &
                                 temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                 grid, rmin, n_irr, eps_)
                         end if
                      end if
                   end do
                end if
             end if
          end do
       end do
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE fco_26

  !!<summary>Finds the symmetry preserving HNFs for the face centered
  !!orthorhombic lattice with determinant n. Assuming the basis of A =
  !![[ 1.  , 1.  , -1.  ], [-1.779796, 0.1798 , 0.  ], [ 0.735376,
  !!-1.61953 , -1.68415 ]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE fco_16(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta11, beta12, gamma11, gamma12, gamma21
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f,c)==0) then
          do e=0, (f-1), c
             gamma21 = -2*e+e*e/c
             if (MOD(gamma21, f)==0) then
                do b=0, (c-1)
                   beta12 = -a+2*b
                   if (MOD(beta12, c)==0) then
                      do d=0, (f-1)
                         if ((MOD(e*beta11, c)==0) .and. (MOD(e*beta12, c)==0)) then
                            beta11 = beta12-d
                            gamma11 = -e*beta11/c
                            gamma12 = 2*d-e*beta12/c
                            if ((MOD(beta11, c)==0) .and. (MOD(gamma11, f)==0) .and. &
                                 (MOD(gamma12, f)==0)) then
                               nhnfs = nhnfs + 1          
                               if (all_hnfs) then 
                                  temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                               else
                                  
                                  temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                                  call compare_grids(lat_vecs, B_vecs, at, &
                                       temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                       grid, rmin, n_irr, eps_)
                               end if
                            end if
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE fco_16
    
  !!<summary>Finds the symmetry preserving HNFs for the body centered
  !!orthorhombic lattice with determinant n. Assuming the basis of A =
  !![[0.5,1,1.5],[0,2,0],[0,0,3]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE bco_19(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma12, gamma13, beta13
    integer :: es(2), nes
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, 2)==0) then
          nes = 2
          es = (/0, (f/2)/)
       else
          nes = 1
          es = (/0, -1/)
       end if

       do b = 0, (c-1)
          beta13 = a +2*b
          if (MOD(beta13, c)==0) then
             do k=1, nes
                e = es(k)
                if (MOD(e*beta13, c)==0) then
                   gamma13 = e*beta13/c
                   if (MOD(gamma13, f)==0) then
                      do d=0, (f-1)
                         gamma12 = a + 2*d
                         if (MOD(gamma12, f)==0) then
                            nhnfs = nhnfs + 1          
                            if (all_hnfs) then 
                               temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                            else
                               
                               temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                               call compare_grids(lat_vecs, B_vecs, at, &
                                    temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                    grid, rmin, n_irr, eps_)
                            end if
                         end if
                      end do
                   end if
                end if
             end do
          end if
       end do
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE bco_19

  !!<summary>Finds the symmetry preserving HNFs for the body centered
  !!orthorhombic lattice with determinant n. Assuming the basis of A =
  !![[ 1.41144 , 0.0885622, -2.  ], [-0.99868 , 2 .21232 , 1.268178 ],
  !![ 3.41012 , -1.1237578, -1.268178 ]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE bco_8(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta21, beta11, gamma21, gamma11, gamma13
    integer :: bs(2), nbs
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(2*f, c)==0) then
          if (MOD(c,2.0_dp)==0) then
             nbs = 2
             bs = (/0, (c/2)/)
          else
             nbs = 1
             bs = (/0, -1/)
          end if
          
          do j=1, nbs
             b = bs(j)
             do e=0, (f-1)
                beta21 = 2*e
                if (MOD(beta21*e, c)==0) then
                   gamma21 = -beta21+beta21*e/c
                   if ((MOD(beta21, c)==0) .and. (MOD(gamma21, f)==0)) then
                      do d=0, (f-1)
                         beta11 = -a+(2*b)-2*d
                         if (MOD(beta11*e, c)==0) then
                            gamma11 = -beta11*e/c
                            gamma13 = a+(2*d)-b*beta21/c
                            if (all_hnfs) then 
                               temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                            else
                                        
                               temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                               call compare_grids(lat_vecs, B_vecs, at, &
                                    temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                    grid, rmin, n_irr, eps_)
                            end if
                         end if
                      end do
                   end if
                end if
             end do
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE bco_8

  !!<summary>Finds the symmetry preserving HNFs for the body centered
  !!orthorhombic lattice with determinant n. Assuming the basis of A =
  !![[-1.53633, 1.36706, -1.33073], [ 1.  , 1.  , 1.  ], [ 1.61803,
  !!-0.61803, -1.  ]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE bco_42(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta11, gamma12, gamma11, gamma13
    integer :: es(2), nes
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, 2)==0) then
          nes = 2
          es = (/0, (f/2)/)
       else
          nes = 1
          es = (/0, -1/)
       end if

       do j=1, nes
          e = es(j)
          do b=0, (c-1)
             beta11 = -a+2*b
             if (MOD(beta11*e, c)==0) then
                gamma12 = -beta11*e/c
                if ((MOD(beta11, c)==0) .and. (MOD(gamma12, f)==0)) then
                   do d=0, (f-1)
                      gamma11 = -a+(2*d)-e*beta11/c
                      gamma13 = -a+2*d
                      if ((MOD(gamma11, f)==0) .and. (MOD(gamma13, f)==0)) then
                         nhnfs = nhnfs + 1          
                         if (all_hnfs) then 
                            temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                         else
                            
                            temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                            call compare_grids(lat_vecs, B_vecs, at, &
                                 temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                 grid, rmin, n_irr, eps_)
                         end if
                      end if
                   end do
                end if
             end if
          end do
       end do

    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE bco_42
       
  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!orthorhombic lattice with determinant n. Assuming A =
  !![[0.5,1,0],[0.5,-1,0],[0,0,3]] (1st basis choince in
  !!notes/base_ortho) for case 38 and A = [[ 1.  , 1.  , 1.  ], [ 1.
  !!, -1.  , -1.  ], [ 0.  , -1.73205, 1.73205]] for case
  !!13.</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE baseco_38_13(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma13, gamma23, beta13
    integer :: es(2), ds(2), nes, nds
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(c,a)==0) then
          if (MOD(f, 2)==0) then
             nes = 2
             nds = 2
             es = (/0, (f/2)/)
             ds = (/0, (f/2)/)
          else
             nes = 1
             nds = 1
             es = (/0, -1/)
             ds = (/0, -1/)
          end if

          b = 0
          do while (b<c)
             if (MOD(b, a)==0) then
                beta13 = -a +b*b/a
                if (MOD(beta13, c)==0) then
                   do j = 1, nes
                      e = es(j)
                      do k = 1, nds
                         d = ds(k)
                         if (MOD(b*d, a)==0) then
                            gamma13 = -d + b*d/a -e*beta13/c
                            gamma23 = c*d/a -e -b*e/a
                            if ((MOD(gamma13, f)==0) .and. (MOD(gamma23, f)==0)) then
                               nhnfs = nhnfs + 1          
                               if (all_hnfs) then 
                                  temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                               else
                                  
                                  temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                                  call compare_grids(lat_vecs, B_vecs, at, &
                                       temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                       grid, rmin, n_irr, eps_)
                               end if
                            end if
                         end if
                      end do
                   end do
                end if
             end if
             b = b + a
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE baseco_38_13

  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!orthorhombic lattice with determinant n. Assuming A =
  !![[-0.3333333, -1.54116 , 1.87449 ], [ 1.  , 1.  , 1.  ], [ 2.  ,
  !!-1.  , -1.  ]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE baseco_23(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta13, gamma13, gamma23
    integer :: es(2), bs(2), nes, nbs
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, a)==0) then
          if (MOD(c, 2)==0) then
             nbs = 2
             bs = (/0, (c/2)/)
          else
             nbs = 1
             bs = (/0, -1/)
          end if
          if (MOD(f, 2)==0) then
             nes = 2
             es = (/0, (f/2)/)
          else
             nes = 1
             es = (/0, -1/)
          end if

          do j=1, nbs
             b = bs(j)
             if (MOD(b*f, a*c)==0) then
                do k=1, nes
                   e = es(k)
                   if ((MOD(b*e, a*c)==0) .and. (MOD(2*b*e, c*f)==0) .and. (MOD(e, a)==0)) then
                      do d=0, (f-1), a
                         if ((MOD(b*d, a)==0)) then
                            beta13 = -b+b*d/a
                            if (MOD(beta13*e, c)==0 .and. (MOD(d*d, a)==0)) then
                               gamma13 = -a+(d*d/a)-beta13*e/c
                               gamma23 = e+(d*e/a)-b*e*e/(a*c)
                               if ((MOD(beta13, c)==0) .and. (MOD(gamma13, f)==0) .and. &
                                    (MOD(gamma23, f)==0)) then
                                  nhnfs = nhnfs + 1          
                                  if (all_hnfs) then 
                                     temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                          0, c, e, &
                                          0, 0, f/),(/3,3/))
                                  else
                                     
                                     temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                          0, c, e, &
                                          0, 0, f/),(/3,3/))
                                     call compare_grids(lat_vecs, B_vecs, at, &
                                          temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                          grid, rmin, n_irr, eps_)
                                  end if
                               end if
                            end if
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE baseco_23

  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!orthorhombic lattice with determinant n. Assuming A = [[ 1.  , 1.
  !!, 1.  ], [ 1.61803 , -0.618034, -1.  ], [-1.05557 , 1.99895 ,
  !!-0.943376]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE baseco_40(n,spHNFs)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta13, gamma13, gamma12, gamma22
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, c)==0) then
          do e=0, (f-1), c
             gamma22 = 2*e-e*e/c
             if (MOD(gamma22, f)==0) then
                do d=0, (f-1), c
                   gamma12 = 2*d-d*e/c
                   if (MOD(gamma12, f)==0) then
                      do b=0, (c-1)
                         beta13 = 2*b-d
                         gamma13 = beta13*e/c
                         if ((MOD(beta13, c)==0) .and. (MOD(gamma13, f)==0)) then
                            nhnfs = nhnfs + 1          
                            if (all_hnfs) then 
                               temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                            else
                               
                               temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                               call compare_grids(lat_vecs, B_vecs, at, &
                                    temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                    grid, rmin, n_irr, eps_)
                            end if
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE baseco_40

  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!orthorhombic lattice with determinant n. Assuming A = [[1, 1, 1],
  !![1.41421, -1.41421, 0], [-1.43541, -1.43541, 1.37083]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE baseco_36(n,spHNFs)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta12, beta22, beta32, gamma13, gamma12, gamma22
    integer :: bs(2), es(2), nes, nbs
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, a)==0) then
          if (MOD(f, 2)==0) then
             nes = 2
             es = (/0, (f/2)/)
          else
             nes = 1
             es = (/0, -1/)
          end if
          if (MOD(c, 2)==0) then
             nbs = 2
             bs = (/0, (c/2)/)
          else
             nbs = 1
             bs = (/0, -1/)
          end if

          do j=1, nbs
             b = bs(j)
             beta32 = -b*f/a
             if (MOD(beta32, c)==0) then
                do k=1, nes
                   e = es(k)
                   if ((MOD(b*e, a)==0) .and. (MOD(b*e, c)==0)) then
                      beta22 = -b*e/a
                      gamma13 = -2*b*e/c
                      if ((MOD(beta22, c)==0) .and. (MOD(gamma13, f)==0) .and. &
                           (MOD(e, a)==0)) then
                         do d=0, (f-1)
                            if ((MOD(b*d, a)==0) .and. (MOD(d*d, a)==0) .and. &
                                 (MOD(d*e, a)==0)) then
                               beta12 = -b*d/a
                               if (MOD(beta12*e, c)==0) then
                                  gamma12 = 2*d-(d*d/a)-beta12*e/c
                                  gamma22 = 2*e-(d*e/a)-beta22*e/c
                                  if ((MOD(beta12, c)==0) .and. (MOD(gamma12, f)==0) &
                                       .and. (MOD(gamma22, f)==0)) then
                                     nhnfs = nhnfs + 1          
                                     if (all_hnfs) then 
                                        temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                     else
                                        
                                        temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                        call compare_grids(lat_vecs, B_vecs, at, &
                                             temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                             grid, rmin, n_irr, eps_)
                                     end if
                                  end if
                               end if
                            end if
                         end do
                      end if
                   end if
                end do
             end if
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE baseco_36
    
  !!<summary>Finds the symmetry preserving HNFs for the simple
  !!monoclinic lattice with determinant n. Assuming the basis of A =
  !![[2,0,0],[0,2,0],[0.5,0,2]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE sm_33(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma12
    integer :: es(2), bs(2), nes, nbs
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(c, 2)==0) then
          nbs = 2
          bs = (/0, (c/2)/)
       else
          nbs = 1
          bs = (/0, -1/)
       end if
       if (MOD(f, 2)==0) then
          nes = 2
          es = (/0, (f/2)/)
       else
          nes = 1
          es = (/0, -1/)
       end if

       do j=1, nbs
          b = bs(j)
          do k=1, nes
             e = es(k)
             if (MOD(2*b*e, c) == 0) then
                gamma12 = 2*b*e/c
                if (MOD(gamma12, f)==0) then
                   do d=0, (f-1)
                      nhnfs = nhnfs + 1          
                      if (all_hnfs) then 
                         temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                      else
                                     
                         temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                         call compare_grids(lat_vecs, B_vecs, at, &
                              temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                              grid, rmin, n_irr, eps_)
                      end if
                   end do
                end if
             end if
          end do
       end do
       deallocate(es,bs)
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE sm_33
    
  !!<summary>Finds the symmetry preserving HNFs for the simple
  !!monoclinic lattice with determinant n. Assuming the basis of A =
  !![[1,1,1],[1.61803,-0.618034,-1],[-0.668912,1.96676,-1.29785]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE sm_35(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma12
    integer :: bs(2), nbs
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(c, 2)==0) then
          nbs = 2
          bs = (/0, (c/2)/)
       else
          nbs = 1
          bs = (/0, -1/)
       end if

       do j=1, nbs
          b = bs(j)
          do d=0, (f-1)
             do e=0, (f-1)
                if (MOD(2*b*e, c)==0) then
                   gamma12 = 2*d-2*b*e/c
                   if (MOD(gamma12, f)==0) then
                      nhnfs = nhnfs + 1          
                      if (all_hnfs) then 
                         temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                      else
                                        
                         temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                         call compare_grids(lat_vecs, B_vecs, at, &
                              temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                              grid, rmin, n_irr, eps_)
                      end if
                   end if
                end if
             end do
          end do
       end do
       deallocate(bs)
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE sm_35
    
  !!<summary>Finds the symmetry preserving HNFs for the simple
  !!monoclinic lattice with determinant n. Assuming the basis of A =
  !![[1,1,1],[1.22474487,-1.22474487,-1],[-0.16598509,-1.64308297,1.80906806]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE sm_34(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: ds(2), es(2), nd_es
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, 2)==0) then
          nd_es = 2
          ds = (/0, (f/2)/)
          es = (/0, (f/2)/)
       else
          nd_es = 1
          ds = (/0, -1/)
          es = (/0, -1/)
       end if

       do j=1, nd_es
          d = ds(j)
          do k=1, nd_es
             e = es(k)
             do b=0, (c-1)
                nhnfs = nhnfs + 1          
                if (all_hnfs) then 
                   temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                        0, c, e, &
                        0, 0, f/),(/3,3/))
                else
                                        
                   temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                        0, c, e, &
                        0, 0, f/),(/3,3/))
                   call compare_grids(lat_vecs, B_vecs, at, &
                        temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                        grid, rmin, n_irr, eps_)
                end if
             end do
          end do
       end do
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE sm_34
    
  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!monoclinic lattice with determinant n. Assuming the basis of A =
  !![[1,1,0],[0,2,0],[0.5,0,2]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE basecm_14(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma12, beta12
    integer :: es(2), nes
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, 2)==0) then
          nes = 2
          es = (/0, (f/2)/)
       else
          nes = 1
          es = (/0, -1/)
       end if

       do b=0, (c-1)
          beta12 = a + 2*b
          if (MOD(beta12, c)==0) then
             do k=1, nes
                e = es(k)
                if (MOD(e*beta12, c) == 0) then
                   gamma12 = e*beta12/c
                   if (MOD(gamma12, f)==0) then
                      do d=0, (f-1)
                         nhnfs = nhnfs + 1          
                         if (all_hnfs) then 
                            temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                         else
                            
                            temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                            call compare_grids(lat_vecs, B_vecs, at, &
                                 temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                 grid, rmin, n_irr, eps_)
                         end if
                      end do
                   end if
                end if
             end do
          end if
       end do
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE basecm_14
    
  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!monoclinic lattice with determinant n. Assuming for basis 10 A =
  !![[-1.46391, 0.  , 1.96391], [ 1.  , 1.  , 1.  ], [ 0.  , 2.  , 0.
  !!]], for basis 17 A = [[-0.05387 , -0.61088 , 2.51474 ], [ 1.  , 1.
  !!, 1.  ], [ 1.809568, -0.15957 , 0.  ]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE basecm_10_17(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma12, gamma22
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       do e=0, (f-1)
          gamma22 = c+2*e
          if (MOD(gamma22, f)==0) then
             do d=0, (f-1)
                do b=0, (c-1)
                   gamma12 = b+2*d
                   if (MOD(gamma12, f)==0) then
                      nhnfs = nhnfs + 1          
                      if (all_hnfs) then 
                         temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                      else
                         
                         temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                         call compare_grids(lat_vecs, B_vecs, at, &
                              temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                              grid, rmin, n_irr, eps_)
                      end if
                   end if
                end do
             end do
          end if
       end do
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE basecm_10_17
    
  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!monoclinic lattice with determinant n. Assuming for basis 20 A =
  !![[ 1.  , 1.  , 1.  ], [ 1.70119 , -1.45119 , 1.  ], [ 0.69779 ,
  !!-1.4322505, 3.23446 ]], for basis 25 A = [[ 1.  , 1.  , 1.  ], [ 1
  !!.45119, -1.70119, -1.  ], [ 0.28878, -3.26895,
  !!0.48018]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE basecm_20_25(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma12, gamma22
    integer :: bs(2), nbs
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(c, 2)==0) then
          nbs = 2
          bs = (/0, (c/2)/)
       else
          nbs = 1
          bs = (/0, -1/)
       end if

       do e=0, (f-1)
          gamma22 = -c-2*e
          if (MOD(gamma22, f)==0) then
             do k=1, nbs
                b = bs(k)
                if (MOD(2*b*e, c)==0) then
                   gamma12 = -b-2*b*e/c
                   if (MOD(gamma12, f)==0) then
                      do d=0, (f-1)
                         nhnfs = nhnfs + 1          
                         if (all_hnfs) then 
                            temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                         else
                            
                            temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                            call compare_grids(lat_vecs, B_vecs, at, &
                                 temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                 grid, rmin, n_irr, eps_)
                         end if
                      end do
                   end if
                end if
             end do
          end if
       end do
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE basecm_20_25
    
  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!monoclinic lattice with determinant n. Assuming A = [[ 0.464824,
  !!-1.464824, -1.907413], [-1.618033, 0.618033, -1.  ], [-1.  , -1.
  !!, 0.  ]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE basecm_27(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma12, gamma22
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       do e=0, (f-1)
          gamma22 = c+2*e
          if (MOD(gamma22, f)==0) then
             do b=0, (c-1)
                do d=0, (f-1)
                   gamma12 = a+b+2*d
                   if (MOD(gamma12, f)==0) then
                      nhnfs = nhnfs + 1          
                      if (all_hnfs) then 
                         temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                      else
                                        
                         temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                         call compare_grids(lat_vecs, B_vecs, at, &
                              temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                              grid, rmin, n_irr, eps_)
                      end if
                   end if
                end do
             end do
          end if
       end do
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE basecm_27
    
  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!monoclinic lattice with determinant n. Assuming A = [[-1.44896 ,
  !!0.948958, -1.  ], [-1.  , -1.  , 0.  ], [ 0.342424, -1.342424,
  !!-2.02006 ]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE basecm_28(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma21, gamma11
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, c)==0) then
          do e=0, (f-1), c
             gamma21 = 2*e+e*e/c
             if (MOD(gamma21, f)==0) then
                do d=0, (f-1), c
                   gamma11 = 2*d+d*e/c
                   if (MOD(gamma11, f)==0) then
                      do b=0, (c-1)
                         nhnfs = nhnfs + 1          
                         if (all_hnfs) then 
                            temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                         else
                                        
                            temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                            call compare_grids(lat_vecs, B_vecs, at, &
                                 temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                 grid, rmin, n_irr, eps_)
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE basecm_28
    
  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!monoclinic lattice with determinant n. Assuming for basis 29 A =
  !![[-0.666125, 1.16613 , 2.04852 ], [ 1.  , 1.  , 0.  ], [ 1.61803 ,
  !!-0.618034, 1.  ]], for basis 30 A = [[ 1.  , 1.  , 0.  ], [
  !!1.61803 , -0.618034 , 1.  ], [-0.0361373, 0.536137 , 2.38982
  !!]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE basecm_29_30(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma21, gamma11
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, c)==0) then
          do e=0, (f-1), c
             gamma21 = 2*e+e*e/c
             if (MOD(gamma21, f)==0) then
                do d=0, (f-1), c
                   gamma11 = 2*d+d*e/c
                   if (MOD(gamma11, f)==0) then
                      do b=0, (c-1)
                         nhnfs = nhnfs + 1          
                         if (all_hnfs) then 
                            temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                         else
                            
                            temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                            call compare_grids(lat_vecs, B_vecs, at, &
                                 temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                 grid, rmin, n_irr, eps_)
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE basecm_29_30
    
  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!monoclinic lattice with determinant n. Assuming A = [[-1.  , 0.  ,
  !!-1.  ], [ 1.85397 , 0.854143, -1.35397 ], [-1.  , 1.41421 , 1.
  !!]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE basecm_41(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: gamma21, gamma11
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       do e=0, (f-1)
          gamma21 = -c+2*e
          if (MOD(gamma21, f)==0) then
             do b=0, (c-1)
                do d=0, (f-1)
                   gamma11 = -b+2*d
                   if (MOD(gamma11, f)==0) then
                      nhnfs = nhnfs + 1          
                      if (all_hnfs) then 
                         temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                      else
                         
                         temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                         call compare_grids(lat_vecs, B_vecs, at, &
                              temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                              grid, rmin, n_irr, eps_)
                      end if
                   end if
                end do
             end do
          end if
       end do
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE basecm_41
    
  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!monoclinic lattice with determinant n. Assuming for basis 37 A =
  !![[-1.79092 , -1.47209 , 0.790922], [ 1.  , 0.  , 1.  ], [ 1.  ,
  !!-1.41421 , -1.  ]], for basis 39 A = [[ 0.  , 1.73205 , 1.  ],
  !![-1.  , 0.  , -1.  ], [ 1.66542 , 0.672857, -1.66542 ]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE basecm_37_39(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta11, gamma11
    integer :: es(2), nes
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, 2)==0) then
          nes = 2
          es = (/0, (f/2)/)
       else
          nes = 1
          es = (/0, -1/)
       end if

       do j=1, nes
          e = es(j)
          do b=0, (c-1)
             beta11 = -a+2*b
             if (MOD(beta11*e, c)==0) then
                gamma11 = -beta11*e/c
                if ((MOD(beta11, c)==0) .and. (MOD(gamma11, f)==0)) then
                   do d=0, (f-1)
                      nhnfs = nhnfs + 1          
                      if (all_hnfs) then 
                         temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                      else
                         
                         temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                              0, c, e, &
                              0, 0, f/),(/3,3/))
                         call compare_grids(lat_vecs, B_vecs, at, &
                              temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                              grid, rmin, n_irr, eps_)
                      end if
                   end do
                end if
             end if
          end do
       end do
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE basecm_37_39
    
  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!monoclinic lattice with determinant n. Assuming A = [[-0.39716,
  !!-0.34718, 2.49434], [ 2.64194, -0.14194, 0.  ], [-1.39716,
  !!-1.34718, 1.49434]].</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE basecm_43(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f
    integer :: nds, i, nhnfs, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:)

    integer :: beta12, gamma12, gamma22
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    call get_HNF_diagonals(n,diagonals)

    nds = size(diagonals,2)
    nhnfs = 0
    total_hnfs = 0
    rmin = 0.0_dp
    n_irr = 0
    grid = 0.0_dp

    if (all_hnfs) then 
       do i = 1,nds
          total_hnfs = total_hnfs + diagonals(2,i)*diagonals(3,i)**2
       end do

       allocate(temp_HNFs(3,3,total_hnfs), STAT=status)
       if (status/=0) stop "Failed to allocate memory in hex_12."
    else
       allocate(temp_HNFs(3,3,1))
    end if
    
    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(f, c)==0) then
          do e=0, (f-1), c
             gamma22 = 2*e-e*e/c
             if (MOD(gamma22, f)==0) then
                do d=0, (f-1)
                   beta12 = a+d
                   gamma12 = 2*a+2*d-beta12*e/c
                   if ((MOD(beta12, c)==0) .and. (MOD(gamma12, f)==0)) then
                      do b=0, (c-1)
                         nhnfs = nhnfs + 1          
                         if (all_hnfs) then 
                            temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                         else
                            
                            temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                            call compare_grids(lat_vecs, B_vecs, at, &
                                 temp_HNFs(:,:,1), No, Nu, Co, Cu, O, &
                                 grid, rmin, n_irr, eps_)
                         end if
                      end do
                   end if
                end do
             end if
          end do
       end if
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else 
       allocate(spHNFs(3,3,1))

       spHNFs(:,:,1) = temp_HNFs(:,:,1)
    end if
    
  end SUBROUTINE basecm_43
       
  !!<summary>Finds the symmetry preserving HNFs for the triclinic
  !!lattice with determinant n. Subroutine taken form enumlib on 7/20/17.</summary>
  !!<parameter name="n" regular="true">The target determinant of the
  !!HNFs.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="U" regular="true">Users basis vectors.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="rmin" regular="true">R_min of best HNF.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  SUBROUTINE tric_31_44(n, No, Nu, Co, Cu, O, U, B_vecs, at, spHNFs, grid, rmin, &
       n_irr, nhnfs, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: rmin, grid(3,3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_
    
    integer, pointer    :: d(:,:) => null()
    integer             :: i, j, k, l    ! Loop counters
    integer             :: Nds, Nhnf, ihnf ! # of triplets, # of HNF matrices, HNF counter
    integer             :: status
    real(dp) :: eps
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if
    
    call get_HNF_diagonals(n,d)
    Nds = size(d,2)
    
    ! Count the total number of HNF matrices for given determinant (n)

    if (all_hnfs) then 
       Nhnf = 0
       do i = 1,Nds
          Nhnf = Nhnf + d(2,i)*d(3,i)**2
       enddo

       allocate(spHNFs(3,3,Nhnf),STAT=status)
       if(status/=0) stop "Failed to allocate memory in tric"
    else
       allocate(spHNFs(3,3,1),STAT=status)
    end if
    
    ihnf = 0
    do i = 1,Nds ! Loop over the permutations of the diagonal elements of the HFNs
       do j = 0,d(2,i)-1  ! Look over possible values of row 2, element 1
          do k = 0,d(3,i)-1  ! Ditto for row 3, element 1
             do l = 0,d(3,i)-1  ! Ditto for row 3, element 2
                ihnf = ihnf+1 ! Count the HNFs and construct the next one
                if (all_hnfs) then 
                   spHNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                        0, c, e, &
                        0, 0, f/),(/3,3/))
                else
                   
                   spHNFs(:,:,1) = reshape((/ a, b, d, &
                        0, c, e, &
                        0, 0, f/),(/3,3/))
                   call compare_grids(lat_vecs, B_vecs, at, &
                        spHNFs(:,:,1), No, Nu, Co, Cu, O, &
                        grid, rmin, n_irr, eps_)
                end if
             enddo
          enddo
       enddo  ! End loops over values for off-diagonal elements
    enddo ! End loop over all unique triplets of target determinant (n)
    
    if (ihnf /= Nhnf) stop "HNF: not all the matrices were generated...(bug!)"
  end SUBROUTINE tric_31_44
  
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
  !!given size. Subroutine taken from enumlib on 7/18/17.</summary>
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
