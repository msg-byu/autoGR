module autoGR

    use input_structure
    use find_kgrids
    use kpointgeneration, only: generateIrredKpointList, mapKptsIntoBZ
    use num_types

    implicit none
    private

    public :: get_k_grid_from_autoGR

contains

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