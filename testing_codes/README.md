# kgrid_gen
Kpoint generation codes

## To download the code use (using ssh so you will need to have a public key setup):

```
git clone --recursive git@github.com:wsmorgan/kgrid_gen.git
```

## How to compile.
First compile symlib:

```
cd symlib/src
make F90=gfortran
cd ../../
```

To make the executable that counts how many HNFs exist per volume use:

```
make F90=gfortran srHNF.x
```

To make the executable that lists the point group for a lattice in lattice coordinates use:

```
make F90=gfortran pg.x
```

To change which lattice you are working with go into `struct_enum.in`
and after the work 'bulk' replace the next 3 lines of numbers with
your lattice vectors (each lattice vectors goes on its own line, i.e.,
the lattice vectors are the rows of the matrix in this file).