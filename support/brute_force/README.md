# Brute Force Codes
Kpoint generation codes

## Dependencies

To use this code you will need a copy of `symlib` in the `brute_force` directory.

```
git clone git@github.com:msg-byu/symlib.git

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