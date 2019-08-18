!!<summary>Generates the symmetry preserving HNFs for the different
!!crystal lattices.</summary>
Module sp_hnfs
  use grid_utils
  use vector_matrix_utilities, only: matrix_inverse
  use num_types

  implicit none
  private
  public sc_3, fcc_1, bcc_5, hex_12, hex_22, rhom_9_24, rhom_4_2, st_11, st_21, &
       bct_6_7_15_18, so_32, baseco_23, baseco_36, baseco_40, baseco_38_13, &
       bco_19, bco_8, bco_42, fco_26, fco_16, sm_33, sm_34_35, basecm_10_14_17_27_37_39_41, &
       basecm_43, basecm_28, basecm_29_30, basecm_20_25, tric_31_44

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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grids
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE sc_3(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer, allocatable :: cand_HNFs(:,:,:,:)
    integer :: a,b,c,d,e,f
    integer :: nds, i
    integer :: temp_HNFs(3,3,1), ngrids(2)
    real(dp) :: eps
    real(dp) :: supercell(3,3)
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-3
    end if

    allocate(cand_grids(3,3,1,2), cand_HNFs(3,3,1,2))
    cand_grids = 0.0_dp
    cand_HNFs = 0
    temp_HNFs = 0

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
       elseif ((a==c) .and. (((f/a)==2) .and. (MOD(f, a)==0))) then
          b = 0
          d = a
          e = a
          temp_HNFs(:,:,1) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((c==f) .and. (((f/a)==2) .and. (MOD(f, a)==0))) then
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
    if (all_hnfs) then       
       spHNFs(:,:,1) = temp_HNFs(:,:,1)
       if (.not. all(spHNFs == 0)) then
          nhnfs = 1
       else
          nhnfs = 0
       end if
    else 
       call transform_supercell(temp_HNFs(:,:,1), No, Nu, Co, Cu, O, supercell)
       call matrix_inverse(transpose(supercell), cand_grids(:,:,1,1))
       cand_HNFs(:,:,1,1) = temp_HNFs(:,:,1)
       nhnfs = 1
       ngrids = (/1,0/)
       call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
            grid, spHNFs, best_offset, n_irr, symm_flag, eps)
    end if

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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grids
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE fcc_1(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer, allocatable :: cand_HNFs(:,:,:,:)
    integer :: a,b,c,d,e,f
    integer :: nds, i
    integer :: temp_HNFs(3,3,1), ngrids(2)
    real(dp) :: eps
    real(dp) :: supercell(3,3)
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if
    
    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-3
    end if

    allocate(cand_grids(3,3,1,2), cand_HNFs(3,3,1,2))
    cand_grids = 0.0_dp
    cand_HNFs = 0
    temp_HNFs = 0

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
       elseif ((c==f) .and. ((((f/a) == 2) .and. (MOD(f, a) ==0)) .or. (((f/a) == 4)) &
            .and. (MOD(f, a)==0))) then
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
    if (all_hnfs) then
       spHNFs(:,:,1) = temp_HNFs(:,:,1)
       if (.not. all(spHNFs == 0)) then
          nhnfs = 1
       else
          nhnfs = 0
       end if
    else 
       call transform_supercell(temp_HNFs(:,:,1), No, Nu, Co, Cu, O, supercell)
       call matrix_inverse(transpose(supercell), cand_grids(:,:,1,1))
       cand_HNFs(:,:,1,1) = temp_HNFs(:,:,1)
       nhnfs = 1
       ngrids = (/1,0/)
       call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
            grid, spHNFs, best_offset, n_irr, symm_flag, eps)
    end if

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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grids
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE bcc_5(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer, allocatable :: cand_HNFs(:,:,:,:)
    integer :: a,b,c,d,e,f
    integer :: nds, i
    integer :: temp_HNFs(3,3,1), ngrids(2)
    real(dp) :: eps
    real(dp) :: supercell(3,3)
    logical :: all_hnfs

    if (present(all_hnfs_)) then
       all_hnfs = all_hnfs_
    else
       all_hnfs = .False.
    end if
    
    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-3
    end if

    allocate(cand_grids(3,3,1,2), cand_HNFs(3,3,1,2))
    cand_grids = 0.0_dp
    cand_HNFs = 0
    temp_HNFs = 0
    
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
       elseif ((a==c) .and. (((f/a)==2) .and. (MOD(f, a)==0))) then
          b = 0
          d = a
          e = a
          temp_HNFs(:,:,1) = reshape((/ a, b, d, &
               0, c, e, &
               0, 0, f/),(/3,3/))
          exit
       elseif ((a==c) .and. (((f/a)==4) .and. (MOD(f, a)==0))) then
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
    allocate(spHNFs(3,3,1))
    if (all_hnfs) then
       spHNFs(:,:,1) = temp_HNFs(:,:,1)
       if (.not. all(spHNFs == 0)) then
          nhnfs = 1
       else
          nhnfs = 0
       end if
    else 
       call transform_supercell(temp_HNFs(:,:,1), No, Nu, Co, Cu, O, supercell)
       call matrix_inverse(transpose(supercell), cand_grids(:,:,1,1))
       cand_HNFs(:,:,1,1) = temp_HNFs(:,:,1)
       nhnfs = 1
       ngrids = (/1,0/)
       call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
            grid, spHNFs, best_offset, n_irr, symm_flag, eps)
    end if

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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE hex_12(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j, nes
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: beta13, beta11, gamma13, gamma11, gamma12, gamma21, gamma22
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
       eps = 1E-3
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
                if ((MOD(beta13, int(c,li))==0) .and. (MOD(beta11, int(c,li))==0)) then
                   if (MOD(f, 2)==0) then
                      es = (/0, (f/2)/)
                      nes = 2
                   else
                      es = (/0, -1/)
                      nes = 1
                   end if
                   do j = 1, nes
                      e = es(j)
                      gamma13 = (a+2*b)*e
                      if (MOD(gamma13, int(c,li))==0) then
                         gamma13 = gamma13/c
                         if (MOD(gamma13, int(f,li))==0) then
                            do d=0, (f-1)
                               if ((MOD(b*d, a)==0) .and. (MOD(e*beta11, int(c,li))==0) .and. &
                                    (MOD(b*e, a)==0) .and. (MOD(c*d, a)==0) .and. &
                                    (MOD(c*d-b*e, a)==0)) then
                                  gamma11 = b*d/a -e*beta11/c
                                  gamma12 = 2*d + b*d/a - e*beta11/c
                                  gamma21 = c*d/a - 2*e - b*e/a
                                  gamma22 = (c*d - b*e)/a
                                  if ((MOD(gamma11,int(f,li))==0) .and. &
                                       (MOD(gamma12,int(f,li))==0) .and. &
                                       (MOD(gamma21,int(f,li))==0) .and. &
                                       (MOD(gamma22,int(f,li)) ==0)) then
                                     nhnfs = nhnfs + 1
                                     if (all_hnfs) then
                                        temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                     else

                                        temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                        call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                             Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                             ngrids, eps)
                                     end if
                                  end if
                               end if
                            end do
                         end if
                      end if
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE hex_22(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: gamma21, gamma22, gamma11, gamma12
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
       eps = 1E-3
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
                if ((MOD(gamma21, int(f,li))==0) .and. (MOD(gamma22, int(f,li))==0)) then
                   do b=0, (c-1)
                      do d=0, (f-1), c
                         if (MOD(d*e, c)==0) then
                            gamma11 = -b+d*2
                            gamma12 = -b+d-(d*e)/c
                            if ((MOD(gamma11, int(f,li))==0) .and. &
                                 (MOD(gamma12, int(f,li))==0)) then
                               nhnfs = nhnfs + 1
                               if (all_hnfs) then
                                  temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                               else

                                  temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                                  call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                       Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                       ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if
  end SUBROUTINE hex_22

  !!<summary>Finds the symmetry preserving HNFs for the rhombohedral
  !!lattice with determinant n. Assuming the basis of A =
  !![[1,2,2],[2,1,2],[4,3,3]]for basis 9, A = [[-0.255922,-1.44338,0.92259],
  !![1.51184,0,-0.845178],[1.255922,1.44338,0.07741]]
  !!for basis 24.</summary>
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE rhom_9_24(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: beta13, beta22, beta12, gamma11, gamma12, gamma21, gamma22
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
       eps = 1E-3
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

    do i =1, nds
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
             beta13 = f+b*f/a
             if (MOD(beta13,int(c,li))==0) then
                e = 0
                do while (e <f)
                   beta22 = e + b*e/a
                   gamma21 = c + 2*e
                   gamma11 = b + 2*b*e/c
                   if ((MOD(beta22, int(c,li))==0) .and. (MOD(gamma21, int(f,li))==0) .and. &
                        (MOD(gamma11, int(f,li))==0)) then
                      d = 0
                      do while (d < f)
                         beta12 = -a + b + d + d*b/a
                         if (MOD(beta12, int(c,li))==0) then
                            gamma12 = -b -d + (d*d/a) - e*beta12/c
                            gamma22 = -c -2*e + (d*e/a) - e*beta22/c
                            if ((MOD(gamma12, int(f,li))==0) .and. &
                                 (MOD(gamma22, int(f,li))==0)) then
                               nhnfs = nhnfs + 1
                               if (all_hnfs) then
                                  temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                               else
                                  
                                  temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                                  call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                       Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                       ngrids, eps)
                               end if
                            end if
                         end if
                         d = d + a
                      end do
                   end if
                   e = e + a
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if

  end SUBROUTINE rhom_9_24

  !!<summary>Finds the symmetry preserving HNFs for the rhombohedral
  !!lattice with determinant n. Assuming A = [[-1.11652,-0.610985,0.616515],
  !![0.0,-1.32288,-0.5],[1.0,1.32288,1.5]]for basis 2 and A =
  !![[-0.548584,0.774292,1.04858],[0.0,-1.32288,0.5],[1.0,1.32288,0.5]]
  !!for basis 4.</summary>
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE rhom_4_2(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: beta32, beta22, beta12, gamma11, gamma12, gamma22, gamma21
    integer :: nbs
    integer :: bs(2)
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
       eps = 1E-3
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

       if (MOD(f, a)==0)then
          if (MOD(c, 2)==0) then
             nbs = 2
             bs = (/0, (c/2)/)
          else
             nbs = 1
             bs = (/0, -1/)
          end if

          do j=1, nbs
             b = bs(j)
             if (MOD(b*f, a)==0) then
                beta32 = -f+b*f/a
                if(MOD(beta32, int(c,li))==0) then
                   do e=0, (f-1)
                      if ((MOD(b*e, a)==0) .and. (MOD(2*b*e, c)==0)) then
                         beta22 = -e+b*e/a
                         gamma11 = b-2*b*e/c
                         gamma21 = c-2*e
                         if((MOD(e, a)==0) .and. (MOD(beta22, int(c,li))==0) .and. &
                              (MOD(gamma11, int(f,li))==0) .and. &
                              (MOD(gamma21, int(f,li))==0))then
                            do d=0, (f-1)
                               if (MOD(b*d, a)==0) then
                                  beta12 = -a+b-d+b*d/a
                                  if ((MOD(e*beta12, int(c,li))==0) .and. (MOD(d*d, a)==0) &
                                       .and. (MOD(d*e, a)==0)) then
                                     gamma12 = (-a+d*d/a)-(e*beta12/c)
                                     gamma22 = (-e+d*e/a)-(e*beta22/c)
                                     if((MOD(d, a)==0) .and. (MOD(beta12, int(c,li))==0) .and. &
                                          (MOD(gamma12, int(f,li))==0) .and. &
                                          (MOD(gamma22, int(f,li))==0))then
                                        nhnfs = nhnfs + 1
                                        if (all_hnfs) then
                                           temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                                0, c, e, &
                                                0, 0, f/),(/3,3/))
                                        else

                                           temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                                0, c, e, &
                                                0, 0, f/),(/3,3/))
                                           call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                                Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                                ngrids, eps)
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
      end if
   end do

   if (all_hnfs) then
      allocate(spHNFs(3,3,nhnfs))

      spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
   else
      allocate(spHNFs(3,3,1))
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if

  end SUBROUTINE rhom_4_2

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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE st_11(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: beta13, gamma13, gamma12, gamma23
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
       eps = 1E-3
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
                if (MOD(beta13, int(c,li))==0) then
                   do k = 1, nes
                      e = es(k)
                      if (MOD(2*b*e, c)==0) then
                         gamma12 = 2*b*e/c
                         if (MOD(gamma12, int(f,li))==0) then
                            do d = 0, (f-1)
                               if ((MOD(d*b, a)==0) .and. (MOD(d*c, a)==0) &
                                    .and. (MOD(e*b, a)==0)) then
                                  gamma13 = -e*beta13/c + d*(b/a-1)
                                  gamma23 = -e*(b/a+1) +d*c/a
                                  if ((MOD(gamma13, int(f,li))==0) .and. &
                                       (MOD(gamma23, int(f,li))==0)) then
                                     nhnfs = nhnfs + 1
                                     if (all_hnfs) then
                                        temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                     else

                                        temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                        call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                             Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                             ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE st_21(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: beta13, gamma13, gamma12, gamma23
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
       eps = 1E-3
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
             nbs = 1
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
                if (MOD(beta13, int(c,li))==0) then
                   do z=1, nes
                      e = es(z)
                      if ((MOD(2*d*e, c)==0) .and. (MOD(e*e, c)==0)) then
                         gamma12 = 2*d-2*d*e/c
                         gamma13 = -b+d-e*beta13/c
                         gamma23 = -c+e*e/c
                         if ((MOD(gamma12, int(f,li))==0) .and. &
                              (MOD(gamma13, int(f,li))==0) .and. &
                              (MOD(gamma23, int(f,li))==0)) then
                            nhnfs = nhnfs + 1
                            if (all_hnfs) then
                               temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                            else

                               temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                               call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                    Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                    ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if

  end SUBROUTINE st_21

  !!<summary>Finds the symmetry preserving HNFs for the body centered
  !!tetragonal lattice with determinant n. Assuming the basis of
  !!A  = [[1.80278,-1.47253,0.762655],[2.80278,0.13535,-0.791285],
  !![0.80278,-0.47253,2.762655]] for 6,
  !!A = [[1.95095, 1.19163, 0.879663],[0.0, 2.60788, 0.44606],
  !![0.95095, -0.41625, 2.433603]] for 7,
  !!A = [[-1.0,-1.0,2.0],[0.0,-2.0,0.0],[-2.0,0.0,0.0]] for 15,
  !!A = [[-2.0,-1.0,1.0],[-3.0,1.0,0.0],[-1.0,-3.0,0.0]] for 18.</summary>
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE bct_6_7_15_18(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: gamma21, gamma13, beta12, gamma12
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
       eps = 1E-3
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

       if (MOD(f, c)==0)then
          if (MOD(f, 2)==0)then
             nes = 2
             es = (/0, (f/2)/)
          else
             nes = 1
             es = (/0, -1/)
          end if
          do j=1, nes
             e = es(j)
             if (MOD(e, c)==0) then
                gamma21 = -c+e*e/c
                if (MOD(gamma21, int(f,li))==0) then
                   do d=0, (f-1)
                      gamma13 = a+2*d
                      if (MOD(gamma13, int(f,li))==0) then
                         do b=0, (c-1)
                            beta12 = b-d
                            if (MOD(beta12, int(c,li))==0) then
                               gamma12 = -b+d-(e*beta12/c)
                               if(MOD(gamma12, int(f,li))==0)then
                                  nhnfs = nhnfs + 1
                                  if (all_hnfs) then
                                     temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                          0, c, e, &
                                          0, 0, f/),(/3,3/))
                                  else

                                     temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                          0, c, e, &
                                          0, 0, f/),(/3,3/))
                                     call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                          Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                          ngrids, eps)
                                  end if
                               end if
                            end if
                         end do
                      end if
                   end do
                end if
             end if
          end do
       endif
    end do

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs))

       spHNFs(:,:,1:nhnfs) = temp_HNFs(:,:,1:nhnfs)
    else
       allocate(spHNFs(3,3,1))
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if

  end SUBROUTINE bct_6_7_15_18

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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE so_32(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j, k, z
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

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
       eps = 1E-3
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

       if (MOD(c, 2)==0) then
          nbs = 2
          bs = (/0, (c/2)/)
       else
          nbs = 1
          bs = (/0, -1/)
       end if
       if (MOD(f, 2)==0) then
          ne_ds = 2
          es = (/0, (f/2)/)
          ds = (/0, (f/2)/)
       else
          ne_ds = 1
          es = (/0, -1/)
          ds = (/0, -1/)
       end if

       do j = 1, nbs
          b = bs(j)
          do k=1, ne_ds
             e = es(k)
             if (MOD((2*b*e), (f*c))==0) then
                do z = 1, ne_ds
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
                      call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                           Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                           ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE fco_26(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: gamma12, gamma13, gamma23
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in fco_26."
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
                if ((MOD(gamma23, int(f,li))==0) .and. (MOD(gamma12, int(f,li))==0)) then
                   do d = 0, (f-1)
                      gamma13 = a +b +2*d
                      if (MOD(gamma13, int(f,li))==0) then
                         nhnfs = nhnfs + 1
                         if (all_hnfs) then
                            temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                         else

                            temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                            call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                 Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                 ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if

  end SUBROUTINE fco_26

  !!<summary>Finds the symmetry preserving HNFs for the face centered
  !!orthorhombic lattice with determinant n. Assuming the basis of A =
  !![[1.04442, 1.43973, 1.68415], [0.779796, -1.1789, 1.0],
  !![1.779796, -0.1789, 0]].</summary>
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE fco_16(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: gamma11, gamma12, gamma21
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in fco_16."
    else
       allocate(temp_HNFs(3,3,1))
    end if

    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(c, 2)==0) then
          nbs = 2
          bs = (/(c/2), 0/)
       else
          nbs = 1
          bs = (/0, -1/)
       end if
       do j=1, nbs
          b = bs(j)
          do e=0, (f-1)
             if (MOD(2*b*e, c)==0) then
                gamma11 = -b-2*b*e/c
                gamma21 = c+2*e
                if (MOD(gamma11, int(f,li))==0 .and. MOD(gamma21, int(f,li))==0) then
                   do d=0, (f-1)
                      gamma12 = a+2*d-(2*b*e/c)
                      if (MOD(gamma12, int(f,li))==0) then
                         nhnfs = nhnfs + 1
                         if (all_hnfs) then
                            temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                         else

                            temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                            call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                 Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                 ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE bco_19(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: gamma12, gamma13, beta13
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in bco_19."
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
          if (MOD(beta13, int(c,li))==0) then
             do k=1, nes
                e = es(k)
                if (MOD(e*beta13, int(c,li))==0) then
                   gamma13 = e*beta13/c
                   if (MOD(gamma13, int(f,li))==0) then
                      do d=0, (f-1)
                         gamma12 = a + 2*d
                         if (MOD(gamma12, int(f,li))==0) then
                            nhnfs = nhnfs + 1
                            if (all_hnfs) then
                               temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                            else

                               temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                               call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                    Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                    ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE bco_8(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: beta21, beta11, gamma21, gamma11, gamma13
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in bco_8."
    else
       allocate(temp_HNFs(3,3,1))
    end if

    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(2*f, c)==0) then
          if (MOD(c, 2)==0) then
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
                if (MOD(beta21*e, int(c,li))==0) then
                   gamma21 = -beta21+beta21*e/c
                   if ((MOD(beta21, int(c,li))==0) .and. (MOD(gamma21, int(f,li))==0)) then
                      do d=0, (f-1)
                         beta11 = -a+(2*b)-2*d
                         if (MOD(beta11, int(c,li))==0) then
                            gamma11 = -beta11*e/c
                            gamma13 = a+(2*d)-b*beta21/c
                            if ((MOD(gamma11, int(f,li))==0) .and. &
                                 (MOD(gamma13, int(f,li))==0)) then
                               nhnfs = nhnfs + 1
                               if (all_hnfs) then
                                  temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                               else

                                  temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                                  call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                       Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                       ngrids, eps)
                               end if
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE bco_42(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: beta11, gamma12, gamma11, gamma13
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in bco_42."
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
             if (MOD(beta11*e, int(c,li))==0) then
                gamma12 = -beta11*e/c
                if ((MOD(beta11, int(c,li))==0) .and. (MOD(gamma12, int(f,li))==0)) then
                   do d=0, (f-1)
                      gamma11 = -a+(2*d)-e*beta11/c
                      gamma13 = -a+2*d
                      if ((MOD(gamma11, int(f,li))==0) .and. (MOD(gamma13, int(f,li))==0)) then
                         nhnfs = nhnfs + 1
                         if (all_hnfs) then
                            temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                         else

                            temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                 0, c, e, &
                                 0, 0, f/),(/3,3/))
                            call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                 Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                 ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE baseco_38_13(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: gamma13, gamma23, beta13
    integer :: es(2), ds(2), ne_ds
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in baseco_38_13."
    else
       allocate(temp_HNFs(3,3,1))
    end if

    do i =1,nds
       a = diagonals(1,i)
       c = diagonals(2,i)
       f = diagonals(3,i)

       if (MOD(c,a)==0) then
          if (MOD(f, 2)==0) then
             ne_ds = 2
             es = (/0, (f/2)/)
             ds = (/0, (f/2)/)
          else
             ne_ds = 1
             es = (/0, -1/)
             ds = (/0, -1/)
          end if

          b = 0
          do while (b<c)
             if (MOD(b*b, a)==0) then
                beta13 = -a +b*b/a
                if (MOD(beta13, int(c,li))==0) then
                   do j = 1, ne_ds
                      e = es(j)
                      do k = 1, ne_ds
                         d = ds(k)
                         if (MOD(b*d, a)==0) then
                            gamma13 = -d + b*d/a -e*beta13/c
                            gamma23 = c*d/a -e -b*e/a
                            if ((MOD(gamma13, int(f,li))==0) .and. &
                                 (MOD(gamma23, int(f,li))==0)) then
                               nhnfs = nhnfs + 1
                               if (all_hnfs) then
                                  temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                               else

                                  temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                       0, c, e, &
                                       0, 0, f/),(/3,3/))
                                  call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                       Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                       ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE baseco_23(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: beta13, gamma13, gamma23
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in baseco_23."
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
                            if (MOD(beta13*e, int(c,li))==0 .and. (MOD(d*d, a)==0)) then
                               gamma13 = -a+(d*d/a)-beta13*e/c
                               gamma23 = e+(d*e/a)-b*e*e/(a*c)
                               if ((MOD(beta13, int(c,li))==0) .and. &
                                    (MOD(gamma13, int(f,li))==0) .and. &
                                    (MOD(gamma23, int(f,li))==0)) then
                                  nhnfs = nhnfs + 1
                                  if (all_hnfs) then
                                     temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                          0, c, e, &
                                          0, 0, f/),(/3,3/))
                                  else

                                     temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                          0, c, e, &
                                          0, 0, f/),(/3,3/))
                                     call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                          Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                          ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE baseco_40(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: beta13, gamma13, gamma12, gamma22
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in baseco_40."
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
             if (MOD(gamma22, int(f,li))==0) then
                do d=0, (f-1), c
                   gamma12 = 2*d-d*e/c
                   if (MOD(gamma12, int(f,li))==0) then
                      do b=0, (c-1)
                         beta13 = 2*b-d
                         gamma13 = beta13*e/c
                         if ((MOD(beta13, int(c,li))==0) .and. &
                              (MOD(gamma13, int(f,li))==0)) then
                            nhnfs = nhnfs + 1
                            if (all_hnfs) then
                               temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                            else

                               temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                    0, c, e, &
                                    0, 0, f/),(/3,3/))
                               call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                    Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                    ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE baseco_36(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: beta12, beta22, beta32, gamma13, gamma12, gamma22
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in baseco_36."
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
             if (MOD(beta32, int(c,li))==0) then
                do k=1, nes
                   e = es(k)
                   !The (MOD(e, a)==0) check is an actual condition
                   !that comes from the integer relationships. Not a
                   !mistake in droping the b from b*e.
                   if ((MOD(e, a)==0) .and. (MOD(b*e, c)==0)) then
                      beta22 = -b*e/a
                      gamma13 = -2*b*e/c
                      if ((MOD(beta22, int(c,li))==0) .and. (MOD(gamma13, int(f,li))==0)) then
                         do d=0, (f-1)
                            if ((MOD(b*d, a)==0) .and. (MOD(d*d, a)==0)) then
                               beta12 = -b*d/a
                               if (MOD(beta12*e, int(c,li))==0) then
                                  gamma12 = 2*d-(d*d/a)-beta12*e/c
                                  gamma22 = 2*e-(d*e/a)-beta22*e/c
                                  if ((MOD(beta12, int(c,li))==0) .and. &
                                       (MOD(gamma12, int(f,li))==0) &
                                       .and. (MOD(gamma22, int(f,li))==0)) then
                                     nhnfs = nhnfs + 1
                                     if (all_hnfs) then
                                        temp_HNFs(:,:,nhnfs) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                     else

                                        temp_HNFs(:,:,1) = reshape((/ a, b, d, &
                                             0, c, e, &
                                             0, 0, f/),(/3,3/))
                                        call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                             Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                             ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE sm_33(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: gamma12
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in sm_33."
    else
       allocate(temp_HNFs(3,3,1))
    end if

    do i =1, nds
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
                if (MOD(gamma12, int(f,li))==0) then
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
                         call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                              Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                              ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if

  end SUBROUTINE sm_33

  !!<summary>Finds the symmetry preserving HNFs for the simple
  !!monoclinic lattice with determinant n. Assuming the basis of A =
  !![[1,1,1],[1.22474487,-1.22474487,-1],[-0.16598509,-1.64308297,1.80906806]]
  !!for 34, and a =  =[[-0.668912,1.96676,-1.29785],[1.61803,-0.618034,-1.0]
  !!,[1.0,1.0,1.0]] for 35.</summary>
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE sm_34_35(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j, k
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in sm_34_35."
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
                   call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                        Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                        ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if

  end SUBROUTINE sm_34_35

  !!<summary>Finds the symmetry preserving HNFs for the base centered
  !!monoclinic lattice with determinant n. Assuming:
  !!for basis 10 A = [[1, -1, 1],[-1.46391, 0, 1.96391],[0, 2, 0]],
  !!for basis 14 A = [[-1,1,0],[0.5,0,2],[0,-2,0]],
  !!for basis 17 A = [[-1.05387,-1.61088,1.51474],[-0.244302,-2.77045,0.51474],[1.809568,-0.15957,0.0]]
  !!for basis 27 A = [[-1.464824,0.464824,1.907413],[-0.153209,0.153209,-2.907413],[1.0,1.0,0.0]],
  !!for basis 37 A = [[-1.79092,-1.47209,0.790922],[1.0,-1.41421,-1.0],[1.0,0.0,1.0]],
  !!for basis 39 A = [[0, -1.73205,-1],[-1.66542, -0.672857, 1.66542], [1,0,1]],
  !!for basis 41 A = [[-1.85397, -0.854143, 1.35397],[1, 0, 1],[1, -1.41421, -1]].</summary>
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE basecm_10_14_17_27_37_39_41(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status, j
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: gamma11
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in basecm_10_14_17_27_37_39_41."
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
          do d=0, (f-1)
             e = es(j)
             gamma11 = -1*a + 2*d
             if (MOD(gamma11,int(f,li))==0) then
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
                      call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                           Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                           ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if

  end SUBROUTINE basecm_10_14_17_27_37_39_41

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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE basecm_20_25(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, k, status
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: gamma12, gamma22
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in basecm_20_25."
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
          if (MOD(gamma22, int(f,li))==0) then
             do k=1, nbs
                b = bs(k)
                if (MOD(2*b*e, c)==0) then
                   gamma12 = -b-2*b*e/c
                   if (MOD(gamma12, int(f,li))==0) then
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
                            call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                 Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                 ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if

  end SUBROUTINE basecm_20_25

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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE basecm_28(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: gamma21, gamma11
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in basecm_28."
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
             if (MOD(gamma21, int(f,li))==0) then
                do d=0, (f-1), c
                   gamma11 = 2*d+d*e/c
                   if (MOD(gamma11, int(f,li))==0) then
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
                            call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                 Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                 ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE basecm_29_30(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: gamma21, gamma11
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in basecm_29_30."
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
             if (MOD(gamma21, int(f,li))==0) then
                do d=0, (f-1), c
                   gamma11 = 2*d+d*e/c
                   if (MOD(gamma11, int(f,li))==0) then
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
                            call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                 Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                 ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if

  end SUBROUTINE basecm_29_30

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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE basecm_43(n, No, Nu, Co, Cu, O, U, B_vecs, at, offsets, best_offset, spHNFs, grid, &
       n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, pointer :: diagonals(:,:) => null()
    integer :: a,b,c,d,e,f, best_HNF(3,3)
    integer :: nds, i, status
    integer(li) :: total_hnfs
    integer, allocatable :: temp_HNFs(:,:,:), cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer(li) :: beta12, gamma12, gamma22
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
       eps = 1E-3
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
       if (status/=0) stop "Failed to allocate memory in basecm_43."
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
             if (MOD(gamma22, int(f,li))==0) then
                do d=0, (f-1)
                   beta12 = a+d
                   gamma12 = 2*a+2*d-beta12*e/c
                   if ((MOD(beta12, int(c,li))==0) .and. (MOD(gamma12, int(f,li))==0)) then
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
                            call compare_grids(temp_HNFs(:,:,1), No, Nu, &
                                 Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                                 ngrids, eps)
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
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
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
  !!<parameter name="B_vecs">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="mult" regular="true">Multiplication factor to
  !!reach correct determinant.</parameter>
  !!<parameter name="offsets" regular="true">The possible offsets to
  !!apply.</parameter>
  !!<parameter name="spHNFs" regular="true">The symmetry preserving
  !!HNFs.</parameter>
  !!<parameter name="grid" regular="true">The kpoint grid
  !!found.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="nhnfs" regular="true">The number of HNFs found.</parameter>
  !!<parameter name="offsets" regular="true">The offsets to try with
  !!the grids.</parameter>
  !!<parameter name="best_offset" regular="true">The offset that works
  !!with the best grid.</parameter>
  !!<parameter name="all_hnfs_" regular="true">True if all HNFs are
  !!wanted.</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  SUBROUTINE tric_31_44(n, No, Nu, Co, Cu, O, U, B_vecs, at, mult, offsets, best_offset, &
       spHNFs, grid, n_irr, nhnfs, symm_flag, eps_, all_hnfs_)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: spHNFs(:,:,:)
    real(dp), allocatable :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: Co(3,3), Cu(3,3), mult, symm_flag
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3), U(3,3), offsets(:,:)
    integer, intent(out) :: n_irr, nhnfs
    real(dp), intent(out) :: grid(3,3), best_offset(3)
    logical, optional, intent(in) :: all_hnfs_
    real(dp), optional, intent(in) :: eps_

    integer, allocatable :: cand_HNFs(:,:,:,:)
    real(dp), allocatable :: cand_grids(:,:,:,:)
    integer :: ngrids(2)
    real(dp) :: rmin(2)

    integer, pointer    :: d(:,:) => null()
    integer             :: i, j, k, l    ! Loop counters
    integer             :: Nds, ihnf ! # of triplets, HNF counter
    integer             :: status
    integer :: best_HNF(3,3)
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
       eps = 1E-3
    end if

    call get_HNF_diagonals(n,d)
    Nds = size(d,2)

    ! Count the total number of HNF matrices for given determinant (n)
    nhnfs = 0
    do i = 1,Nds
      nhnfs = nhnfs + d(2,i)*d(3,i)**2
    enddo

    if (all_hnfs) then
       allocate(spHNFs(3,3,nhnfs),STAT=status)
       if(status/=0) stop "Failed to allocate memory in tric"
    else
       allocate(spHNFs(3,3,1),STAT=status)
    end if

    rmin = 0.0_dp
    ihnf = 0
    do i = 1,Nds ! Loop over the permutations of the diagonal elements of the HFNs
       do j = 0,d(2,i)-1  ! Look over possible values of row 2, element 1
          do k = 0,d(3,i)-1  ! Ditto for row 3, element 1
             do l = 0,d(3,i)-1  ! Ditto for row 3, element 2
                ihnf = ihnf+1 ! Count the HNFs and construct the next one
                if (all_hnfs) then
                   spHNFs(:,:,ihnf) = reshape((/ d(1,i), j, k, &
                        0, d(2,i), l, &
                        0, 0, d(3,i)/),(/3,3/))
                   spHNFs(:,:,ihnf) = spHNFs(:,:,ihnf)*mult
                else

                   spHNFs(:,:,1) = reshape((/ d(1,i), j, k, &
                        0, d(2,i), l, &
                        0, 0, d(3,i)/),(/3,3/))
                   spHNFs(:,:,1) = spHNFs(:,:,1)*mult
                   call compare_grids(spHNFs(:,:,1), No, Nu, &
                        Co, Cu, O, cand_grids, rmin, cand_HNFs, &
                        ngrids, eps)
                end if
             enddo
          enddo
       enddo  ! End loops over values for off-diagonal elements
    enddo ! End loop over all unique triplets of target determinant (n)

    if (ihnf /= nhnfs) stop "HNF: not all the matrices were generated...(bug!)"
    if (.not. all_hnfs) then
       if (any(ngrids > 0)) then
          call grid_selection(U, B_vecs, at, cand_grids, cand_HNFs, ngrids, offsets, &
               grid, best_HNF, best_offset, n_irr, symm_flag, eps)
          if (any(best_HNF > 0)) then
             spHNFs(:,:,1) = best_HNF
          else
             spHNFs = 0
          end if
       else
          spHNFs = 0
       end if
    end if
  end SUBROUTINE tric_31_44

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
