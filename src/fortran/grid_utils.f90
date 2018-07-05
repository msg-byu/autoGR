Module grid_utils
  use num_types
  use vector_matrix_utilities, only: matrix_inverse, minkowski_reduce_basis, norm, cross_product, determinant
  use numerical_utilities, only: equal
  use kpointgeneration, only: generateIrredKpointList

  implicit none
  private
  public transform_supercell, grid_selection, grid_metrics, compare_grids

CONTAINS

  !!<summary>Finds the symmetry preserving supercell in the users
  !!basis.</summary>
  !!<parameter name="spHNF" regular="true">One of the symmetry
  !!preserving HNFs in our basis</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="UB" regular="true">The supercell in the users
  !!basis.</parameter>
  SUBROUTINE transform_supercell(spHNF,No,Nu,Co,Cu,O,UB)
    integer, intent(in) :: spHNF(3,3)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3)
    real(dp), intent(out) :: UB(3,3)

    integer :: F(3,3)
    real(dp) :: Oinv(3,3), Noinv(3,3), Cuinv(3,3), L(3,3)

    call matrix_inverse(O,Oinv)
    call matrix_inverse(No,Noinv)
    call matrix_inverse(real(Cu,dp),Cuinv)
    L = matmul(matmul(O,spHNF),Oinv)
    F = nint(matmul(Noinv,matmul(matmul(L,O),Co)))
    UB = matmul(matmul(Nu,F),Cuinv)

  end SUBROUTINE transform_supercell

  !!<summary>Finds the rmin and the number of irreducible k-points for
  !!lattice given an HNF.</summary>
  !!<parameter name="lat_vecs" regular="true">The lattice vectors for
  !!the crystal.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="HNF" regular="true">The list of generating
  !!vectors for the candidate grids.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="rmin" regular="true">The rmin for this grid.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducbile
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="grid" regular="true">The grid from the HNF.</parameter>
  SUBROUTINE grid_metrics(lat_vecs, B_vecs, at, HNF, No, Nu, Co, Cu, O, grid, rmin, n_irr, eps_)
    real(dp), intent(in) :: lat_vecs(3,3)
    real(dp), optional, intent(in) :: eps_
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: HNF(3,3)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3)
    real(dp), intent(out) :: rmin, grid(3,3)
    integer, intent(out) :: n_irr

    real(dp) :: supercell(3,3), shift(3), norms(3), reduced_grid(3,3)
    real(dp) :: grid_inv(3,3), lat_trans(3,3)
    real(dp) :: eps
    real(dp)              :: R(3,3)
    real(dp), pointer     :: rdKlist(:,:)
    integer, pointer      :: weights(:)

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    shift = 0.0_dp

    call transform_supercell(HNF, No, Nu, Co, Cu, O, supercell)

    grid_inv = transpose(supercell)
    call matrix_inverse(grid_inv, grid)
    lat_trans = transpose(lat_vecs)
    call matrix_inverse(lat_trans, R)

    call minkowski_reduce_basis(grid, reduced_grid, eps)
    norms(1) = sqrt(dot_product(reduced_grid(:,1), reduced_grid(:,1)))
    norms(2) = sqrt(dot_product(reduced_grid(:,2), reduced_grid(:,2)))
    norms(3) = sqrt(dot_product(reduced_grid(:,3), reduced_grid(:,3)))
    rmin = min(norms(1), norms(2), norms(3))
    call generateIrredKpointList(lat_vecs, B_vecs, at, grid, R, shift, rdKlist, weights, eps_=eps)
    n_irr = size(rdKlist,1)

  end SUBROUTINE grid_metrics

  !!<summary>Compares a grid to the metrics of another to see which is
  !!has the better rmin or number of irreducible k-points.</summary>
  !!<summary>Finds the rmin and the number of irreducible k-points for
  !!lattice given an HNF.</summary>
  !!<parameter name="lat_vecs" regular="true">The lattice vectors for
  !!the crystal.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="HNF" regular="true">The list of generating
  !!vectors for the candidate grids.</parameter>
  !!<parameter name="No" regular="true">Our niggli basis.</parameter>
  !!<parameter name="Nu" regular="true">The users niggli
  !!basis.</parameter>
  !!<parameter name="Co" regular="true">Transformation from our basis
  !!to the our niggli basis.</parameter>
  !!<parameter name="Cu" regular="true">Transformation from the users
  !!basis to the users niggli basis.</parameter>
  !!<parameter name="O" regular="true">Our basis vectors.</parameter>
  !!<parameter name="rmin" regular="true">The rmin for these grids.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="grids" regular="true">A list of grids with the best
  !!r_min.</parameter>
  !!<parameter name="best_hnfs" regular="true">The HNFs for the grids with
  !!this r_min.</parameter>
  !!<parameter name="ngrids" regular="true">The number of grids
  !!currently stored in the list of grids.</parameter>
  SUBROUTINE compare_grids(lat_vecs, B_vecs, at, HNF, No, Nu, Co, Cu, O, grids, rmin, best_HNFs, ngrids, eps_)
    real(dp), intent(in) :: lat_vecs(3,3)
    real(dp), optional, intent(in) :: eps_
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: HNF(3,3)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3)
    real(dp), intent(inout) :: rmin(2)
    real(dp), allocatable, intent(inout) :: grids(:,:,:,:)
    integer, intent(inout) :: ngrids(2)
    integer, allocatable, intent(inout) :: best_HNFs(:,:,:,:)

    real(dp) :: supercell(3,3), shift(3), reduced_grid(3,3), norms(3), supercell2(3,3), red_supercell(3,3)
    real(dp) :: lat_trans(3,3)
    real(dp) :: eps
    real(dp)              :: R(3,3)
    integer, allocatable :: ralloc_HNFs(:,:,:,:)
    real(dp), allocatable :: ralloc_grids(:,:,:,:)
    real(dp) :: temp_rmin, temp_grid(3,3), temp_grid_inv(3,3)

    integer :: j

    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-3
    end if

    shift = 0.0_dp
    if (.not. allocated(grids)) then
       allocate(grids(3,3,10,2))
       grids = 0.0_dp
       ngrids = 0
    end if
    if (.not. allocated(best_HNFs)) then
       allocate(best_HNFs(3,3,10,2))
       best_HNFs = 0
    end if

    call transform_supercell(HNF, No, Nu, Co, Cu, O, supercell)

    temp_grid_inv = transpose(supercell)
    call matrix_inverse(transpose(supercell), temp_grid)
    lat_trans = transpose(lat_vecs)
    call matrix_inverse(lat_trans, R)

    call minkowski_reduce_basis(temp_grid, reduced_grid, eps)
    norms(1) = sqrt(dot_product(reduced_grid(:,1), reduced_grid(:,1)))
    norms(2) = sqrt(dot_product(reduced_grid(:,2), reduced_grid(:,2)))
    norms(3) = sqrt(dot_product(reduced_grid(:,3), reduced_grid(:,3)))
    temp_rmin = min(norms(1), norms(2), norms(3))
    if (temp_rmin > (rmin(1)+rmin(1)*eps) .and. (.not. equal(temp_rmin, rmin(1), eps))) then
       rmin(2) = rmin(1)
       ngrids(2) = ngrids(1)
       grids(:,:,:,2) = 0.0_dp
       best_HNFs(:,:,:,2) = 0
       best_HNFs(1:3,1:3,1:ngrids(2),2) = best_HNFs(1:3,1:3,1:ngrids(2),1)
       grids(1:3,1:3,1:ngrids(2),2) = grids(1:3,1:3,1:ngrids(2),1)
       
       rmin(1) = temp_rmin
       ngrids(1) = 1
       grids(:,:,:,1) = 0.0_dp
       best_HNFs(:,:,:,1) = 0
       grids(:,:,ngrids(1),1) = temp_grid
       best_HNFs(:,:,ngrids(1),1) = HNF

    else if ((temp_rmin > (rmin(2)+rmin(2)*eps)) .and. (.not. equal(temp_rmin, rmin(1), eps)) &
         .and. (.not. equal(temp_rmin, rmin(2), eps))) then
       rmin(2) = temp_rmin
       ngrids(2) = 1
       grids(:,:,:,2) = 0.0_dp
       best_HNFs(:,:,:,2) = 0
       grids(:,:,ngrids(2),2) = temp_grid
       best_HNFs(:,:,ngrids(2),2) = HNF
    else if (equal(temp_rmin, rmin(1), eps)) then
       if (size(grids,3) == ngrids(1)) then
          allocate(ralloc_HNFs(3,3,ngrids(1),2), ralloc_grids(3,3,ngrids(1),2))
          ralloc_HNFs = best_HNFs
          ralloc_grids = grids
          deallocate(best_HNFs, grids)
          allocate(best_HNFs(3,3,2*ngrids(1),2), grids(3,3,2*ngrids(1),2))
          best_HNFs = 0
          grids = 0.0_dp
          best_HNFs(1:3,1:3,1:ngrids(1),1) = ralloc_HNFs(1:3,1:3,1:ngrids(1),1)
          grids(1:3,1:3,1:ngrids(1),1) = ralloc_grids(1:3,1:3,1:ngrids(1),1)
          best_HNFs(1:3,1:3,1:ngrids(2),2) = ralloc_HNFs(1:3,1:3,1:ngrids(2),2)
          grids(1:3,1:3,1:ngrids(2),2) = ralloc_grids(1:3,1:3,1:ngrids(2),2)
       end if
       ngrids(1) = ngrids(1) + 1
       grids(:,:,ngrids(1),1) = temp_grid
       best_HNFs(:,:,ngrids(1),1) = HNF
    else if (equal(temp_rmin, rmin(2), eps)) then
       if (size(grids,3) == ngrids(2)) then
          allocate(ralloc_HNFs(3,3,ngrids(2),2), ralloc_grids(3,3,ngrids(2),2))
          ralloc_HNFs = best_HNFs
          ralloc_grids = grids
          deallocate(best_HNFs, grids)
          allocate(best_HNFs(3,3,2*ngrids(2),2), grids(3,3,2*ngrids(2),2))
          best_HNFs = 0
          grids = 0.0_dp
          best_HNFs(1:3,1:3,1:ngrids(1),1) = ralloc_HNFs(1:3,1:3,1:ngrids(1),1)
          grids(1:3,1:3,1:ngrids(1),1) = ralloc_grids(1:3,1:3,1:ngrids(1),1)
          best_HNFs(1:3,1:3,1:ngrids(2),2) = ralloc_HNFs(1:3,1:3,1:ngrids(2),2)
          grids(1:3,1:3,1:ngrids(2),2) = ralloc_grids(1:3,1:3,1:ngrids(2),2)
       end if
       ngrids(2) = ngrids(2) + 1
       grids(:,:,ngrids(2),2) = temp_grid
       best_HNFs(:,:,ngrids(2),2) = HNF
    end if

  end SUBROUTINE compare_grids

  !!<summary>Selects the best grid from the list of grids.</summary>
  !!<parameter name="lat_vecs" regular="true">The lattice vectors for
  !!the crystal.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="cand_grids" regular="true">The list of generating
  !!vectors for the candidate grids.</parameter>
  !!<parameter name="best_grid" regular="true">The best grid given the
  !!criteria.</parameter>
  !!<parameter name="ngrids" regular="true">The number of
  !!grids.</parameter>
  !!<parameter name="cand_HNFs" regular="true">The cand HNFs.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducible kpoints
  !!of best_grid</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  SUBROUTINE grid_selection(lat_vecs, B_vecs, at, cand_grids, cand_HNFs, ngrids, best_grid, &
    best_HNF, n_irr, eps_)
    real(dp), intent(in) :: lat_vecs(3,3)
    real(dp), allocatable, intent(in) :: cand_grids(:,:,:,:)
    real(dp), optional, intent(in) :: eps_
    real(dp), intent(out) :: best_grid(3,3)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in), allocatable :: cand_HNFs(:,:,:,:)
    integer, intent(out) :: best_HNF(3,3),n_irr
    integer, intent(in) :: ngrids(2)

    real(dp) :: temp_grid(3,3), norms(3)
    integer :: i, n_irreducible(sum(ngrids)), n_ir_min(1), count
    real(dp) :: eps, det

    real(dp)              :: R(3,3), invLat(3,3), shift(3)
    real(dp), pointer     :: rdKlist(:,:)
    integer, pointer      :: weights(:)
    
    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    shift = 0.0_dp

    call matrix_inverse(lat_vecs, invLat)
    R = transpose(invLat)

    n_irreducible = 0
    count = 0
    do i=1,ngrids(1)
       temp_grid = cand_grids(:,:,i,1)
       det = determinant(temp_grid)
       if (((det<eps) .and. (det >1E-10)) .or. (equal(det,eps,eps))) then
          call generateIrredKpointList(lat_vecs, B_vecs, at, temp_grid, R, shift, rdKlist, weights, eps_=(eps**2))
       else 
          call generateIrredKpointList(lat_vecs, B_vecs, at, temp_grid, R, shift, rdKlist, weights, eps_=eps)
       end if
       n_irreducible(i) = size(rdKlist,1)
       count = count + 1
    end do

    do i=1,ngrids(2)
       temp_grid = cand_grids(:,:,i,2)
       det = determinant(temp_grid)
       if (((det<eps) .and. (det >1E-10)) .or. (equal(det,eps,eps))) then
          call generateIrredKpointList(lat_vecs, B_vecs, at, temp_grid, R, shift, rdKlist, weights, eps_=(eps**2))
       else 
          call generateIrredKpointList(lat_vecs, B_vecs, at, temp_grid, R, shift, rdKlist, weights, eps_=eps)
       end if
       n_irreducible(count+i) = size(rdKlist,1)
    end do
    
    n_ir_min = minloc(n_irreducible)
    n_irr = n_irreducible(n_ir_min(1))
    if (n_ir_min(1) <= count) then
       best_grid = cand_grids(:,:, n_ir_min(1),1)
       best_HNF = cand_HNFs(:,:, n_ir_min(1),1)
    else
       best_grid = cand_grids(:,:, n_ir_min(1)-count,2)
       best_HNF = cand_HNFs(:,:, n_ir_min(1)-count,2)
    end if

  end SUBROUTINE grid_selection


end Module grid_utils
