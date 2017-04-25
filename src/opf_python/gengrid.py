"""The driver for the opf kpoint selection.

This subroutin requires a POSCAR to be present and a GRIDGEN file to
be present.

The POSCAR should be in the standard VASP format. The GRIDGEN file
should specify the target k-point density.
"""

import os
import numpy as np

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


