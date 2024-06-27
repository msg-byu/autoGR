module autoGR

    use input_structure
    use find_kgrids
    use kpointgeneration, only: generateIrredKpointList, mapKptsIntoBZ
    use num_types

    implicit none
    private

    public :: get_k_grid_from_autoGR

contains

  !!<summary>High-level wrapper around autoGR to use when build 
  !!as library.</summary>
  !!<parameter name="n_atoms_in" regular="true">The number of atoms in the input 
  !!structure.</parameter>
  !!<parameter name="lattice_in" regular="true">The real space lattice
  !!vectors.</parameter>
  !!<parameter name="atom_type_in" regular="true">The type of each atom
  !!in the basis.</parameter>
  !!<parameter name="atom_base_in" regular="true">The atomic basis
  !!vectors.</parameter>
  !!<parameter name="nkpts" regular="true">The number of k-points to
  !!search for.</parameter>
  !!<parameter name="offset" regular="true">The offset for the
  !!k-points grid.</parameter>
  !!<parameter name="find_offset" regular="true">'True' if the offset
  !!needs to be determined by the algorithm.</parameter>
  !!<parameter name="best_grid" regular="true">The best k-point grid
  !!found.</parameter>
  !!<parameter name="best_offset" regular="true">The best offset from
  !!those checked for the final grid.</parameter>
  !!<parameter name="lat_id" regular="true">The ID of the niggli reduced
  !!bravais lattice.</parameter>
  !!<parameter regular="true" name="IrrKpList"> List of symmetry-reduced k-points in
  !!Cartesian coordinates. </parameter>
  !!<parameter regular="true" name="weights"> "Weights" of k-points (length of each orbit).
  !!</parameter>
  !!<parameter name="symm_flag" regular="true">Flag that indicates the
  !!symmetries to use.</parameter>
  !!<parameter name="min_kpts_flag" regular="true">Flag that indicates
  !!that the grid with the minimum number of k-points should be
  !!selected rather than the grid with the best folding
  !!ratio.</parameter>
  subroutine get_k_grid_from_autoGR(n_atoms_in,  lattice_in, atom_type_in, &
        atom_base_in, nkpts, offset, find_offset, best_grid, best_offset, &
        lat_id, IRKps, weights, symm_flag, min_kpts)

    integer, intent(in) :: n_atoms_in
    real(dp), intent(in) :: lattice_in(3,3)
    integer,  intent(in) :: atom_type_in(n_atoms_in)
    real(dp), intent(in) :: atom_base_in(3,n_atoms_in)
    real(dp), intent(in) :: offset(3)
    logical, intent(in) :: find_offset
    integer, optional, intent(in) :: symm_flag
    logical, optional, intent(in) :: min_kpts

    integer, intent(inout) :: nkpts

    real(dp), intent(out) :: best_grid(3,3), best_offset(3)
    integer, intent(out) :: lat_id

    real(dp), pointer, intent(out) :: IRKps(:,:)
    integer, pointer, intent(out) :: weights(:)
    
    call  init_input_structure(n_atoms_in,  lattice_in, atom_type_in, atom_base_in)
    
    call find_grid(nkpts, offset, find_offset, best_grid, best_offset, &
        lat_id, symm_flag_=symm_flag, min_kpts_=min_kpts, reps_=reps, aeps_=aeps)

    call generateIrredKpointList(lattice, atom_base, atom_type, best_grid, reduced_R, best_offset, &
        IRKps, weights, reps_=reps, aeps_=aeps, symm_=symm_flag)

    call mapKptsIntoBZ(r_vecs, IRKps, reps)
    
    nkpts = size(IRKps,1)
    
  end subroutine get_k_grid_from_autoGR
    
end module autoGR