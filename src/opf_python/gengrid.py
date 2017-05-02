"""The driver for the opf kpoint selection.

This subroutin requires a POSCAR to be present and a GRIDGEN file to
be present.

The POSCAR should be in the standard VASP format. The GRIDGEN file
should specify the target k-point density.
"""

import os
import numpy as np
from universal import find_srBs

if not os.path.isfile("POSCAR"):
    print("ERROR: The POSCAR file must be provided.")
    exit()

if not os.path.isfile("GRIDGEN"):
    print("ERROR: The GRIDGEN file must be provided.")
    exit()

#Find A from the POSCAR
with open("POSCAR","r") as f:
    f.readline()
    f.readline()

    for i in range(3):
        A.append(float(i) for i in  f.readline().strip().split())

    A = np.transpose(A)

#Read the k-point density from file.
with open("GRIDGEN","r") as f:
    target = int(f.readline().strip())

Bs = find_srBs(A,target)

#now we want to find the grid vectors
grids = []
ns = []
for B in Bs:
    grids.append(np.transpose(np.linalg.inv(B)))
    ns.apppend(np.linalg.det(B)/np.linalg.det(A))
    
#now I want to find the packing fraction and the r_min for each of these grids.
from phenum.vector_utils import _minkowski_reduce_basis

rmin = []
pack_frac = []

for grid in grids:
    min_grid = _minkowski_reduce_basis(grid)
    rm = min(np.linalg.norm(min_grid,axis=0))
    pf = 4/3.*np.pi*rm**2/np.dot(min_grid[0],np.cross(min_grid[1],min_grid[2]))
    rmin.append(rm)
    pack_frac.append(pf)

