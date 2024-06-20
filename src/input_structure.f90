module input_structure
    use num_types
    use vector_matrix_utilities, only: matrix_inverse, determinant, minkowski_reduce_basis
    implicit none
    private

    public :: init_input_structure

    integer, public :: n_atoms
    real(dp), public :: lattice(3,3), r_vecs(3,3), reduced_R(3,3)
    integer, allocatable, public :: atom_type(:)
    real(dp), allocatable, public :: atom_base(:,:)
    real(dp), public :: reps = 1E-6_dp
    real(dp), public :: aeps = 5E-4_dp
    real(dp), public :: r_vol

contains

  subroutine init_input_structure(n_atoms_in,  lattice_in, atom_type_in, atom_base_in)

    integer, intent(in) :: n_atoms_in
    real(dp), intent(in) :: lattice_in(3,3)
    integer,  intent(in) :: atom_type_in(n_atoms_in)
    real(dp), intent(in) :: atom_base_in(3,n_atoms_in)
    real(dp) :: pi= 3.1415926535897932385_dp

    n_atoms = n_atoms_in
    lattice = lattice_in

    allocate(atom_type(n_atoms))
    allocate(atom_base(3,n_atoms))

    atom_type = atom_type_in
    atom_base = atom_base_in

    call matrix_inverse(transpose(lattice),r_vecs)
    call minkowski_reduce_basis(r_vecs, reduced_R, reps)

    r_vol = abs(determinant(2*pi*r_vecs))
    
  end subroutine init_input_structure
    
end module input_structure