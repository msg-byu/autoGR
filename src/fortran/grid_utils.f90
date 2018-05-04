Module grid_utils
  use num_types
  use vector_matrix_utilities, only: matrix_inverse, minkowski_reduce_basis, norm, cross_product
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
  !!<parameter name="rmin" regular="true">The rmin for this grid.</parameter>
  !!<parameter name="n_irr" regular="true">The number of irreducbile
  !!k-points.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  !!<parameter name="grid" regular="true">The grid being compared to
  !!is input. The best grid is output.</parameter>
  SUBROUTINE compare_grids(lat_vecs, B_vecs, at, HNF, No, Nu, Co, Cu, O, grid, rmin, n_irr, eps_)
    real(dp), intent(in) :: lat_vecs(3,3)
    real(dp), optional, intent(in) :: eps_
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)
    integer, intent(in) :: HNF(3,3)
    integer, intent(in) :: Co(3,3), Cu(3,3)
    real(dp), intent(in) :: No(3,3), Nu(3,3), O(3,3)
    real(dp), intent(inout) :: rmin, grid(3,3)
    integer, intent(inout) :: n_irr

    real(dp) :: supercell(3,3), shift(3), reduced_grid(3,3), norms(3)
    real(dp) :: lat_trans(3,3)
    real(dp) :: eps
    real(dp)              :: R(3,3)
    real(dp), pointer     :: rdKlist(:,:)
    integer, pointer      :: weights(:)
    real(dp) :: temp_rmin, temp_grid(3,3), temp_grid_inv(3,3)
    integer :: temp_n_irr
    
    if (present(eps_)) then
       eps = eps_
    else
       eps = 1E-6
    end if

    shift = 0.0_dp    
    
    call transform_supercell(HNF, No, Nu, Co, Cu, O, supercell)

    temp_grid_inv = transpose(supercell)
    call matrix_inverse(temp_grid_inv, temp_grid)
    lat_trans = transpose(lat_vecs)
    call matrix_inverse(lat_trans, R)
    
    call minkowski_reduce_basis(temp_grid, reduced_grid, eps)
    norms(1) = sqrt(dot_product(reduced_grid(:,1), reduced_grid(:,1)))
    norms(2) = sqrt(dot_product(reduced_grid(:,2), reduced_grid(:,2)))
    norms(3) = sqrt(dot_product(reduced_grid(:,3), reduced_grid(:,3)))
    temp_rmin = min(norms(1), norms(2), norms(3))
    
    if (temp_rmin > rmin) then
       call generateIrredKpointList(lat_vecs, B_vecs, at, temp_grid, R, shift, rdKlist, weights, eps_=eps)
       n_irr = size(rdKlist,1)
       rmin = temp_rmin
       grid = temp_grid
       
    else if (equal(temp_rmin, rmin, eps)) then
       call generateIrredKpointList(lat_vecs, B_vecs, at, temp_grid, R, shift, rdKlist, weights, eps_=eps)
       temp_n_irr = size(rdKlist,1)
       
       if (temp_n_irr < n_irr) then 
          rmin = temp_rmin
          n_irr = temp_n_irr
          grid = temp_grid          
       end if       
    end if

  end SUBROUTINE compare_grids

  !!<summary>Selects the best grid from the list of grids.</summary>
  !!<parameter name="lat_vecs" regular="true">The lattice vectors for
  !!the crystal.</parameter>
  !!<parameter name="B_vecs" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="at" regular="true">The atom types of each atom in
  !!the basis.</parameter>
  !!<parameter name="grids" regular="true">The list of generating
  !!vectors for the candidate grids.</parameter>
  !!<parameter name="best_grid" regular="true">The best grid given the
  !!criteria.</parameter>
  !!<parameter name="eps_" regular="true">Floating point
  !!tolerance.</parameter>
  SUBROUTINE grid_selection(lat_vecs,B_vecs,at, grids, best_grid, eps_)
    real(dp), intent(in) :: lat_vecs(3,3)
    real(dp), allocatable, intent(in) :: grids(:,:,:)
    real(dp), optional, intent(in) :: eps_
    real(dp), intent(out) :: best_grid(3,3)
    real(dp), pointer :: B_vecs(:,:)
    integer, intent(inout) :: at(:)

    real(dp) :: reduced_grid(3,3), norms(3)
    integer :: i, n_irreducible
    real(dp) :: r_min, r_min_best, eps

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
    
    r_min_best = 0
    n_irreducible = 0
    do i=1,size(grids,3)
       call minkowski_reduce_basis(grids(:,:,i),reduced_grid,eps)
       norms(1) = sqrt(dot_product(reduced_grid(:,1),reduced_grid(:,1)))
       norms(2) = sqrt(dot_product(reduced_grid(:,2),reduced_grid(:,2)))
       norms(3) = sqrt(dot_product(reduced_grid(:,3),reduced_grid(:,3)))
       r_min = min(norms(1),norms(2),norms(3))
       if ((r_min_best==0) .or. (r_min>r_min_best)) then
          r_min_best = r_min
          best_grid = grids(:,:,i)
          call generateIrredKpointList(lat_vecs,B_vecs,at,best_grid, R, shift, rdKlist, weights, eps_=eps)
          n_irreducible = size(rdKlist,1)
       else if (abs(r_min-r_min_best)<eps) then
          call generateIrredKpointList(lat_vecs,B_vecs,at,best_grid, R, shift, rdKlist, weights, eps_=eps)
          if (size(rdKlist,1) < n_irreducible) then
             r_min_best = r_min
             best_grid = grids(:,:,i)
             n_irreducible = size(rdKlist,1)
          end if
       end if
    end do

  end SUBROUTINE grid_selection
  

end Module grid_utils
