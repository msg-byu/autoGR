# Generate Generalized Regular k-point Grids

Implementation of an algorithm to generate _generalized regular_ k-point grids. The algorithm searches over many GR grids that meet a specified k-point density and returns the one with the highest high symmetry reduction and a good packing fraction.
More information can be found in the following references:

* [GR on-the-fly](https://arxiv.org/abs/1902.03257) paper
* [GR vs MP performance](https://www.sciencedirect.com/science/article/pii/S0927025618304105?via%3Dihub) paper
* [K-point folding](https://arxiv.org/abs/1809.10261) paper
* [Optimal Meshes](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.45.13891) by Moreno and Soler
* [K-point Server](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.155109) by Wisesa, McGill and Mueller.

The algorithm works for all crystal classes but it can occasionally fail when the user-provided structural information is "sloppy" resulting in finite precision errors. An update that will fix this known issue is being developed.

## Compiling the code

To compile the executable, first clone this repository:

```
git clone --recursive git@github.com:msg-byu/GRkgridgen.git
```

The `--recursive` is required so that the dependent repos are also cloned.

To compile the executable, navigate to the `src` directory:

```
cd GRkgridgen/src
```

Compile the executable using:

```
make F90=gfortran kpoints.x
```

The makefile is set up so that the dependencies are also compiled.

## Using the code

The code reads two files, 1) a file the contains the lattice vectors
and atomic basis vectors for the crystal (currently only the VASP
`POSCAR` is supported) and 2) a `KPGEN` file that contains the
desired grid density. The code
outputs a `KPOINTS` file containing the irreducible k-points and their
weights. The `KPGEN` file requires one of the following keyword inputs
to be present to determine the number of points that will be needed in
the grid:

- `NKPTS`: the total number of desired k-points.
- `KPDENSITY`: the desired target k-point density.
- `KSPACING`: the linear space desired between k-points (same as CASTEP inputs).
- `KPPRA`: the desired number of k-points per reciprocal atom.

(We recommend `KPDENSITY` or `KSPACING`.)
Only one of these keywords needs to be specified for the code to run.

If a specific offset for the k-point grid is desired then the
following keyword can be used. If this is not specified, the code will select the best shift. If an unshifted grid is desired, specify a zero shift, `SHIFT= 0.0 0.0 0.0`.

- `SHIFT`: the shift away from the origin for the k-point grid
  (expressed in fractions of the reciprocal lattice vectors).


### Example input files

To generate a grid with a k-point density of 5000 k-points/Å^3 and have the code determine the best offset,
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

Once `POSCAR` and `KPGEN` file have been provided by the user, the
`KPOINTS` file can be generated using the `kpoints.x`. The executable
does not have any optional arguments since all relevant options have
been included in the `KPGEN` file.


## Plots

All the data and scripts used to generate the plots and analysis for
the * [GR on-the-fly](https://arxiv.org/abs/1902.03257) paper can be
found in the `paper` folder.
