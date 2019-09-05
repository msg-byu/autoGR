# Revision History

## Revision 0.7.5
- Modified code so that when grids are rudeced by the minkowski method
  a smaller absolute tolerance is used.

## Revision 0.7.4 (WSM)
- Reverted lots of changes to the initial 0.7.0 release.
- Disabled delinter for now.
- Moved offset selection back to original location in `find_grids.f90`.
- Reverted how offsets are selected.
- Reverted behavior for cases with fewer number of k-points so that they
  are not handled differently, this was a mistake.

## Revision 0.7.3 (WSM)
- Changed search domains to adapt with grid size.
- Changed offests to just be permutations of 0.5 and 0.

## Revision 0.7.2 (WSM)
- Moved the location of code that selects possible offsets.
- Changed the code so that it searches for offsets based off the grid
  generating vectors and not the parent lattic.
- Moved space group generating to earlier in the code to reduce number
  of calls to the subroutine.

## Revision 0.7.1 (WSM)
- Added a check for new POSCAR formats.
- Fixed un-initialized variables in delinter.
- Fixed delinter POSCAR writing formats.
- Fixed allocation of ratio array for fining grids with the fewest
  k-points.

## Revision 0.7.0 (WSM)
- Added `RMIN` as an input option to determine the number of k-points.
- Added symmetry options to the inputs that restrict the number of
  symmetry operations used.
- Added an option to force the minimum number of k-points be returned.
- Added code so that if fewer than 500 k-points are wanted then the
  code will loop over all the grids possible and select the one with the
  fewest points.

## Revision 0.6.8 (GLWH)
- Updated `README.md`

## Revision 0.6.7 (WSM)
- Modified the mods in `sphnfs` to compare long ints to long ints so
  the unit tests will run.
- Fixed logic bug in the `tric` `sphnfs` that the unit tests found.
- Fixed some unit tests.

## Revision 0.6.6 (WSM)
- Changed the workflow to catch errors thrown by `kgridGen`.
- Increased the number of cubic cases checked so that if errors occure
  with `kgridGen` we have more than 3 cubic cases to use (this will
  not happen very often but it's there just in case).
- Changed all cases in `sp_hnfs.f90` so that if `grid_selection`
  doesn't return a grid due to errors in `kgridGen` then the code
  handles it properly.
- Added a check at the end of `find_grid` so that if no grids are
  found it prints an error and stops instead of segfaulting.
- Fixed the offsets so that they no longer trigger warnings in
  `kgridGen` and fixed the triclinic cells offsets (it was allocating
  8 instead of 2).

## Revision 0.6.5 (WSM)
- Changed the code workflow so that the packing fraction is used to
  rule out some grids.
- Changed the code workflow so that at least 5 grids are compared for
  non-cubic and non-triclinic cells.
- Cleaned up the code so that it compiles cleanly.

## Revision 0.6.4 (WMS)
- Added a function tha checks that the point group matches the niggli
  id and stops the code if it doesn't.

## Revision 0.6.3 (WSM)
- Fixed most of the finite precision problems we've been having with
  the interface with `kgridGen`.

## Revision 0.6.2 (WSM)
- Fixed multiple allocations that were causing out of bounds errors.
- Reduced the search space of offsets by taking into account the order
  of the vectors that the user provided.
- Started adding OpenMP functionality.
- Fixed calculation of reciprocal lattice volume in control file.

## Revision 0.6.1 (WSM)
- Changed the search intervals for each crystal class.

## Revision 0.6.0 (WSM)
- Changed the contral file name to KPGEN.
- Changed control file to use key words instead of line numbers for
  greater flexibility.
- Added the ability in the code to search over multiple k-point
  offsets.
- Made the users specified offset get passed all the way through the
  code to ensure we get the best folding for that offset.
- Added the offsets to search over by grid type.

## Revision 0.5.5 (WSM)
- Changed the code so that comparisons to the old rmin and new are
  more flexible.

## Revision 0.5.4 (WSM)
- Changed the code so that it will consider the two largest rmin
  values found in order to explore a larger regoin of the search
  space, this is neccessary because in many cases we weren't finding
  the "optimal" grid because it had a slightly smaller rmin but much
  better folding.
- Changed the default eps value to 1E-3 instead of 1E-6, this makes
  more sence is practice and should be sufficient for our purposes.
- Increased the eps value passed to the k-point folding code in cases
  where the determinant of the new grid is approximately the size of
  the eps value, or when it is smaller than eps but still larger than
  1E-10.

## Revision 0.5.3 (WSM)
- Implemented a faster selection implementation in which a minimal set
  of grids is folded.

## Revision 0.5.2 (WSM)
- Implemented changes that John suggested to fix/improve the
  algorithms.
- Debugged new triclinic case.
- Fixed bugs in sm_34_35 inroduced by refactor.
- Fixed bugs in tric_31_44 introduced by refactor.
- Fixed bugs in bco_8 introduced by refactor.
- Fixed bugs in bcc_5 introduced by refactor.
- Fixed bugs in stet_21 introduced by refactor.
- Fixed bugs in fco_16 introduced by refactor.
- Fixed basis for sm_34 in python and fortran niggli id routines.
- Added the python script used to make test output for fortran to the
  repo in testing_codes/make_output.py.
- Added plotting notebook and some plots.
- Reformatted niggli.f90 to make it a little easier to read.
- Fixed a bug in fing_kgrids.f90 in which the count variable wasn't
  getting set.

## Revision 0.5.1 (WSM)
- Fixed lots of compiler errors in sp_hnfs.f90 and grid_utils.f90.
- Added grid_utils.f90 to the Makefile.
- Removed `smallest_prime` from sp_hnfs.f90.
- Refactored find_grids.f90 to match new code format.
- Rewrote driver.
- Fixed a bug in the cubic k-point number finding algorithm that was
  causing it to skip values at large N.
- Made cubic algorithms integer only, somehow they got missed before.
- Fixed allocations and other problems.
- Added the reciprocal lattice calculation to the comparison routines.

## Revision 0.5.0 (WSM)
- Changed sp_hnfs.f90 to use integers instead of reals.
- Changed sp_hnfs.f90 to return only one grid per determinant size
  unless the `all` bool is set to true in which case it returns all
  the HNFs and no grids. Basically the first round of grid selection
  is now done as the HNFs are being constructed.

## Revision 0.4.1 (WSM)
- Fixed a bug in the fortran driver where a negative determinant was
  being used to find the number of k-points wanted.
- Update basis 2 and 4 with better algorithm.

## Revision 0.4.0 (WSM)
- Fixed the driver to work with the new code setup.
- Fixed the issues in the makefile.
- Cleaned up code so it compiles without warnings.

## Revision 0.3.15 (WSM)
- Fixed a bug in the python code where the users basis and our basis
  got switched.
- Fixed integer rounding errors in fortran code.
- Swapped order of d and e placement in the HNF of baseco_38_13 to
  match the python code.
- Fixed base centered monoclinic 28 for fortran code in sphnfs.f90. It
  was looping over the values of b but not placing them in the HNF.
- Fixed a bug in find_grids.f90 that was preventing the code from
  finding HNFs for basis 29.
- Finished implementing code for grid selection.
- Changed the niggli code so that the path taken to reduce the cell is
  not printed with extra white space.

## Revision 0.3.14 (WSM)
-Changed the type of the L matrix in transform_supercell in find_kgrid.f90.
-Implemented unit tests of transform_supercell.

## Revision 0.3.13 (WSM)
- Re-wrote the get_kpd_cubic subroutine in find_kgrds.f90 to match the
  python functions.
- Fixed bugs in find_kgrids.f90 that prevented it from compiling.
- Added subroutine to find_kgrids.f90 find the kpd and ratio to use
  for triclinic cases.
- Fixed bug with floating point division in universal.py.

## Revision 0.3.12 (WSM)
- Removed redundandt call to base_mono_27 in universal.py.
- Fixed basis determination for triclinic cases in universal.py.

## Revision 0.3.11 (WSM)
- Fixed a sign error in the lattice definition for face centered
  orthorhombic cell 16.

## Revision 0.3.10 (WSM)
- Changed how find_volumes works in universal.py so that all valume
  factors are greater than or equal to the desired k-point
  density. Also made it so that in the cumbic cases only 3 values are
  returned and for all other cases 11 values are returned.
- Updated our basis for face centered orthorhombic case 16 in
  niggli_lat_id.py to match the new one we found.
- Changed find_supercells in universal.py so that it only returns
  supercells for the first 5 non-trivial determinant sizes tried.

## Revision 0.3.9 (WSM)
- Fixed some bugs in the find_volumes subroutine in universal.py.
- Added unit tests for some of the functions in universal.py.
- Changed universal.py so that it only keeps the symmetry presereving
  HNFs of a given size if there are more than 1 of them.
- Fixed bug in the volume finding algorithm.
- Removed deprecated code from hex.py.
- Started adding tests of the find_suprcells python algorithm.

## Revision 0.3.8 (WSM)
- Removed the transpose statements from our basis assignments in
  id_cell.
- Replaced eps*D with eps for the tolerance for niggli case 40 in
  niggli_id.

## Revision 0.3.8 (WSM)
- Fixed sign error in get_params subroutine.
- Fixed unit tests outputs for get_params subroutine (python hadn't
  saved enough digits of precision).
- Fixed the assignment of eps in reduce_cell to be related to abs(vol)
  instead of (vol).
- Fixed logical error when checking to see if the reduced cell
  satisfies the niggli conditions in reduce_cell.
- Changed name of niggli_check to condition_check.

## Revision 0.3.7 (WSM)
- Changed M to temp_M in niggli.f90 to fix compiler errors.
- Removed the manually entered parameters in favor of using num_types.
- Fixed various array declaration typos.
- Fixed several syntax errors.
- Changed a, b, and c in id_cell to temp_a, temp_b, temp_c to fix
  compiler errors.
- Fixed several type declarations.
- Removed ```if ((.not. equal(size(IN,1),3,0))
  .or. (.not. equal(size(IN,2),3,0))) stop "Input matrix must be
  3x3."``` from niggli.f90 since it is redudnant with the array
  declaration.
- Changed others


## Revision 0.3.6 (WSW)
- Added the initial unit testse for the fortran niggli reduction.

## Revision 0.3.5 (WSM)
- Fixed bug in bct_15 and added 2 more tests.
- Fixed some minor bugs in niggli.f90.
- Added niggli.xml.

## Revision 0.3.4 (WSM)
- Fixed and debuged basis bco_42, baseco_23, bct_18, basecm_41,
  basecm_20_25, rhom_24, bco_8, st_21, fco_16, sm_34, basecm_29_30,
  basecm_10_17, and basecm_27.
- Fixed typo in documentation for rhom_24 in rhom.py.
- Added test output and input for smallest_prime function.

## Revision 0.3.3 (WSM)
- Fixed and debugged basis rhom_4_2, baseco_40, basecm_43,
  basecm_37_39, bct_6, hex_22, basecm_28, baseco_36, and bct_7.
- Added testing script for easier fortran to python comparison.

## Revision 0.3.2 (WSM)
- Added unit tests files on each of the niggli cell spHNF supercell
  finding subroutines.
- Fixed some initial bugs in spHNFs.f90 so that tests would compile
  and run.

## Revision 0.3.1 (WSM)
- Fixed a bug in the get_kpd routine, with the change to the new
  niggli approach the lattice id for the fcc, bcc, and sc cases
  changed, that change was not reflected in the get_kpd routine.

## Revision 0.3.0 (WSM)
- Added the algorithms for each niggli basis to sp_hnfs.f90.
- Added niggli.f90 to the code (reduces a cell to its niggli form and
  identifies which niggli cell it is along with 'our' basis choince).
- Modified the find_kgrids.f90 routine to use the new niggli basis
  approach.

## Revision 0.2.13 (WSM)
- Added the sm basis to niggli_id.
- Improved the reliability and accuracy of the niggli reduction
  code. It still needs a bit of work on finding an adaptive floating
  point tolerance.

## Revision 0.2.12 (PH)
- Added the simple monoclinic algorithms.

## Revision 0.2.11 (WSM)
- Added correct basis to the output of the niggli lattice
  identification subroutine (niggli_id).
- Fixed typos.
- Added new subroutine to universal.py to find the transformed
  supercells.
- Rewrote and added to basis_check.py to include all the basis and the
  new transformation.
- Fixed a bug in niggli_lat_id.py where the cell was considered
  positive when the angles were effectively but not actually zero. The
  problem was fixed by replace angle >0 by angle-eps>0.

## Revision 0.2.10 (WSM)
- Implemented unit tests for all niggli cases so far designed.
- Renamed trig.py to rhom.py.
- Removed body_tet_srHNFs_2 from body_tet.py as it is no longer useful.
- Removed base_ortho_23_2 from base_ortho.py as it is no longer useful.
- Removed body_ortho_8_2, body_ortho_42_2, and body_ortho_19_2 from
  body_ortho.py as they are no longer useful.
- Removed face_ortho_16_2 from face_ortho.py as it is no longer useful.
- Renamed sc_srHNFs to sc_3, bcc_srHNFs to bcc_5, fcc_srHNFs to fcc_1,
  trig_srHNFs to rhom_9, stet_srHNFs to stet_11, hex_srHNFs to hex_12,
  body_tet_srHNFs to body_tet_15, base_ortho_srHNFs to
  base_ortho_38_13, body_ortho_srHNFs to body_ortho_19,
  face_ortho_srHNFs to face_ortho_26, and base_mono_srHNFs to
  base_mono_14 to match niggli cell numbers and new naming convention.
- Moved test_basis.py to basis_check.py.
- Added .travis.yml to repo.
- Removed lat_type.py from repo since it's been replaced by
  niggli_lat_id.py.
- Implemented unit tests on niggli_lat_id.py.
- Implemented missing niggli basis 36.
- Started implementing new symmetry preserving grid algorithm.

## Revision 0.2.9 (WSM)
- Implemented all the base centered monoclinic cells.

## Revision 0.2.8 (WSM)
- Implemented all basis options for base centered orthorhombic cells.
- Improved/expanded the basis choices for face centered orthorohmbic cells.
- Added option to print the 'path' of the Niggli Reduction.

## Revision 0.2.7 (WSM)
- Fixed an incerroctly labeled niggli cell hexagnoal cell to the correct family.
- Added option to pass G vector directly into niggli_id.
- Added additional basis options to the body centered orthorohmbic cases.
- Added first additional basis for the face centered orthorhombic cases.
- Identified a niggli cell for each of the remaining base centered
  orthorhombic, simple monoclinic, and base centered monoclinic cells
  in the mathematica notebooks.

## Revision 0.2.6 (WSM)
- Fixed an error in the niggli lat id subroutine in which the wrong
  logic was instituted on some angle comparisons.
- Increased the number of iterations allowed by the niggli reduction
  to allow for larger cells.
- Made the default floating point tolerance dependent on cell volume
  in the niggli reduction process.
- Implemented face centered orthorhombic cell number 26.


## Revision 0.2.5 (WSM)
- Fixed cell number 19 though it's still really slow.
- Implemented cell number 42, the last of the body centered
  orthorhombic cells.

## Revision 0.2.4 (WSM)
- Added body centered orthorhrombic cells for niggli cell 8 and 19. 19
  still has bugs in it though and both are really slow.

## Revision 0.2.3 (WSM)
-Added 4 body centered tetragonal cells to the repository.

## Revision 0.2.2 (WSM)
- Fixed bugs in rhombohedral 24 and rhombohedral 9, the former works
  correctly now and the latter had an incoorect basis, it now runs on
  about the same time scale as rhombehdral 2 and 4.
- Added contitions to rhombohedral_2_4 to make it faster.
- Added function to find the symmetry preserving HNFs for both of the
  hexagonal niggli cells.
- Added the second simple tetragonal niggli basis.

## Revision 0.2.1 (WSM)
- Added tolerances to all comparisons in the niggli identification
  script.
- Added tolerances to all comparsions made in the niggli condition
  checking code.
- Added functions to find the symmetry preserving HNFs for all 4
  niggli rhombohedral cells (rhombohedral 24 is wrong though and
  rhombohedral 9 is painfully slow, rhombohedral 2_4 is about 20 times
  slower than the standard trigonal approach).


## Revision 0.2.0 (WSM)
- Added python niggli reduction code and identification script (still
  needs unit testing).

## Revision 0.1.6 (WSM)
-Added second algorithm for generating the body centered tetragonal
 HNFs.
-Added new output to lat_type that may help the basis
 transformation. It now also outputs the order of the basis vectors
 according to size.
-Fixed bug reported in [issue #4](https://github.com/msg-byu/opf_kgrids/issues/4).
-Added a new text file where I'm writing the comparison between methods section.
-Fixed and updated a number of unit tests.

## Revision 0.1.5 (JJ)
- Renamed unit test test_jeremy.py test_sym_grids.py.
- Fixed remaining import calls.

## Revision 0.1.4 (JJ)
- Fixed import calls in a few modules but some are still missing.
- Renamed jeremy.py sym_grids.py.
- Added latex code for the paper we will submit in the directory named
  paper. Completed rough draft of introduction.
- Added notebooks which will test our symmetry preserving grids on
  empirical pseudopotentials in python (must have BZI to run tests).
- Started tests on empirical pseudopotentials in a notebook called
  Search for Figure of Merit.

## Revision 0.1.3
- Fixed a bug in the python code concerning commencerate lattices.

## Revision 0.1.2
- Updated the driver to print the k-points in reciprocal lattice
  coordinates instead of cartesian. Still needs to be tested to see if
  VASP reads it in properly.

## Revision 0.1.1
- Fixed the first few bugs in the make KPOINTS file routine.

## Revision 0.1.0
- Made a driver that makes a kpoints file from the POSCAR and a file
  containing the desired k-point density and offset.
- Implemented the rest of the code needed to find the optimal k_grid.

## Revision 0.0.14
- Added remaining crystal classes to sp_hnf.f90.
- Added unit tests for all crystal classes.

## Revision 0.0.13
- Added hex, trigonal, simple tetragonal and body centered tetragonal
  cases to sp_hnf.f90.

## 0.0.12
- Began implementing the fortran version of the code (created sc, bcc,
  and fcc subroutines).

## 0.0.11
- Implemented faster algorithms for bcc, fcc, sc, and body centered tetragonal cells.
- Fixed integer division in universal.py for python 3.
- Fixed a bug in lat_type.py that was caused by the appearance of
  negative angles from the dot products.
- Fixed jeremey.py. The loop over the grids was indented to far.

## 0.0.10
- Fixed bugs in jeremy.py and in universal.py.
- Added some new unit tests to check for the bugs.
- Added additional supporting info.

## 0.0.9
- Created driver for grid generation src/opf_python/gengrid.py. The
  driver is still incomplete.
- Created new subroutines in universal.py. The first is to find the
  symmetry preserving grids given a lattice and target density. The
  second returns the volume factors that will be used for generating
  the srHNFs. Both subroutines still need work.
- Added number and time scaling plots for face centered cubic and face
  centered orthorhombic.
- Added a module for Jeremy's test to use, jeremy.py.

## 0.0.8
- Fixed some bugs in lat_type.py.
- Impltemented unit tests for lat_type.py

## 0.0.7
- Added python modules for each of the remaining lattice types to the
src/opf_python directory, as well as their unit tests (all modules
agree to with the brute force method for volumes factors of 1 to 500
in the tests except monoclinic lattices which have only been tested on
volumes of 1 to 150 for time reasons).
-Fixed error in notes/body_centered_ortho.nb.
-Added lat_type.py to repo. It's a module that identifies the lattice type.

## 0.0.6
- Added src/python/trig.py to the repo. It finds the srHNFs for the trigonal case.
- Fixed a logic error in the simple, body centered, and face centered
  cubic cases and added the math for the mathematica notebook in the notes for the
  each case.
- Added src/python/stet.py for the simple tetrogonal case as well as a
  mathematica notebook in the notes folder.
- Added a mathematica notebooks for all remaining lattices.
- Added the paper `Space Group Subgroups generated by Sublattice
  Relations` to the notes folder.
- Fixed an infinite loop in trig.py.
- Implemented unit tests on all existing python modules.

## 0.0.5
- Added the 'canonical' parent lattice to the documentation for each crystal class.
- Fixed bugs in the conditions for bcc and fcc. They now find HNFs at
  the correct volumes in the range (1-100).
- Added jupyter-notebook that shows how to identify the lattice type for the code.
- Added plots of the HNF and srHNF scaling the brute force way in
  addition to the notebook that generated the plots.

## 0.0.4
-Renamed the fcc and bcc subroutines so that they are called fcc and
 bcc. Also changed the comments to reflect that change.

## 0.0.3
- Added python code that handles the sc, bcc, and fcc lattices using
  Rod's divisor method.
- Also added an jupyter notebook that shows the timing for the
  generation of the grids for a single determinant at a time. It
  scales really nicely.

## 0.0.2
- Added the testing_codes dir. This directory contains the brute force
  fortran code that is being used to test how large of an volume
  factor we can generate the srHNFs in a given time is and to
  determine which volumes matter mast.
- Added the src directory which will contain code designed for future
  distribution.
- Added the nots directory in which we can store notes from
  collaborations. For now I'm including Rod's divisors writeup and the
  history of or dialog with Dr. Campbell.


## 0.0.1
- Initial commit
