# Generate Generalized Regular k-point Grids

Implementation of an algorithm to generate generalized regular grids
that have a high symmetry reduction and good packing fraction.  It
works for all crystals that are on an ideal lattice but will
occassionally fail on relaxed systems due to rounding errors. An
update that will fix this known issue is being developed.

## Compiling the code

To compile the executable first clone this repository:

```
git clone --recursive git@github.com:msg-byu/GRkgridgen.git
```

This will clone this repository and its dependencies.

To compile the executable navigate to the `src` directory:

```
cd GRkgridgen/src
```

The executable can then be compiled using:

```
make F90=gfortran kpoints.x
```

## Using the code

The code reads two files, 1) a file the contains the lattice vectors
and atomic basis vectors for the crystal (currently only the VASP
`POSCAR` is supported) and 2) a `KPGEN` a file that contains the
desired grid density, in order to generate the k-point grid. The code
outputs a `KPOINTS` file containing the irreducible k-points and their
weights. The `KPGEN` file requires one of the following keyword inputs
to be present to determine the number of points that will be needed in
the grid:

- `NKPTS`: the total number of desired k-points.
- `KPDENSITY`: the desired target k-point density.
- `KSPACING`: the linear space desired between k-points (same as CASTEP inputs).
- `KPPRA`: the desired number of k-points per reciprocal atom.

Only one of these keywords needs to be specified for the code to run.
If a specific offset for the k-point grid is desired then the
following keyword can be used:

- `SHIFT`: the shift away from the origin for the k-point grid
  (expressed in fractions of the reciprocal lattice vectors).

### Example input files

To generate a grid with a k-point density of 5000 k-points per
reciprocal Angstrom cubed and have the code determine the best offset
the `KPGEN` input file would be:

```
KPDENSITY=5000
```

To generate a grid with an offset of [0.5, 0.5, 0.5], expressed as
fractions of the reciprocal lattice vectors, and 500 k-points the
`KPGEN` input file would be:

```
NKPTS=500
SHIFT= 0.5 0.5 0.5
```

### Generating the k-point grid

Once a `POSCAR` and `KPGEN` file have been created in a folder the
`KPOINTS` file can be generated using the `kpoints.x`. The executable
doesn't have any optional arguments since all relevant options have
been included in the `KPGEN` file.