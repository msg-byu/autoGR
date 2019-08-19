# Generate Generalized Regular k-point Grids

Implementation of an algorithm to generate _generalized regular_ k-point grids. The algorithm searches over many GR grids that meet a specified k-point density and returns the one with the highest high symmetry reduction and a good packing fraction.
More information can be found in the following references:

* [Generalized Regular Grids On-The-Fly](https://arxiv.org/abs/1902.03257) describes the algorithm and method of the code in this repository. It describes how the combinatorial explosion of searching over all possible grids is tamed by generating only symmetry-preserving grids.
* [Performance of Generalized Regular Grids vs Monkhorst-Pack Grids](https://www.sciencedirect.com/science/article/pii/S0927025618304105?via%3Dihub) using over 7,000 calculations of different unit cells and different metallice elements that GR grids are 60% more efficient than MP grids on average. Supplements the work of  Wisesa, McGill and Mueller (2016) listed below.
* [Efficient Algorithm For K-Point Folding](https://arxiv.org/abs/1809.10261) discusses in detail an _O(N)_ algorithm for k-point symmetry reduction. This algorithm is essential to enable the on-the-fly method in the paper above.
* [Optimal Meshes](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.45.13891) by Moreno and Soler (1992). The original proposal to use Generalized Regular grids to accelerate DFT calculations.
* [K-point Server](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.155109) by Wisesa, McGill and Mueller (2016). The first practical implementation of Moreno and Soler's idea.

The algorithm works for all crystal classes, but it can fail occasionally when the user-provided structural information is "sloppy", resulting in large finite precision errors. An update that will fix this known issue is being developed.

## Compiling the code

To compile the executable, first clone this repository:

```
git clone --recursive git@github.com:msg-byu/GRkgridgen.git
```

(If you have trouble with the recursive cloning, it is probably because you don't have ssh keys set up on github for your local machine. Try this: [Adding a new SSH key to your GitHub account](https://help.github.com/en/articles/adding-a-new-ssh-key-to-your-github-account).)

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
- `RMIN`: the minimal distance between points in real space desired.

(We recommend `KPDENSITY` or `KSPACING`.)
Only one of these keywords needs to be specified for the code to run.

If a specific offset for the k-point grid is desired then the
following optional keyword can be used. If this is not specified, the code will select the best shift. If an unshifted grid is desired, specify a zero shift, `SHIFT= 0.0 0.0 0.0`.

- `SHIFT`: the shift away from the origin for the k-point grid
  (expressed in fractions of the reciprocal lattice vectors).

The code's default behavior is to find multiple candidate grids then select the grid that has the best ratio of irreducible points to reducible points. Alternatively the following flag can be used to have the code return the grid that has the fewest number of k-points:

```
MIN_IRR_KPTS = TRUE
```

- `MIN_IRR_KPTS`: when set to `TRUE` the grid with the fewest irreducible k-points will be returned.

It is also possible to restrict the symmetry operations used to fold the k-points using the `USE_SYMMETRY' key word.

- `USE_SYMMETRY`: can take the key words `ALL`, `STRUCTURAL`, `TIME_REVERSAL`, `NONE` where each key word indicates the level of symmetry to restrict the symmetry operations used to fold the k-point grid. `All`, the default, indicates that all the symmetry operations of the reciprocal cell should be used, `STRUCTURAL` indicates that time reversal symmetry should be excluded, `TIME_REVERSAL` indicates that the spacial symmetry operations should be ommited and only time reversal symmetry should be used, `NONE` indicates that the grid should not be folded at all.

### Example input files

To generate a grid with a k-point density of 5000 k-points/Å<sup>^-3</sup> and have the code determine the best offset,
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
the [Generalized Regular Grids On-The-Fly](https://arxiv.org/abs/1902.03257) paper can be
found in the `paper` folder.
