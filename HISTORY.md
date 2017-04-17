# Revision History

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