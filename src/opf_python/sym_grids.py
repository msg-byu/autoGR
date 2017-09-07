"""Finds the possible symmetry preserving grids for a given lattice.
"""

import numpy as np
from opf_python.universal import find_srBs
from opf_python.lat_type import lat_type

def spGrids(A,kpd):
    """Finds the symmetry preserving grids for a given lattice at a given density.

    Args:
        lattice (list of list): The parent lattice vectors as columns of a matrix.
        kpd (int): The target k-point density.
    
    Returns:
        grids (list of dict): A list of dictionaries contaning the grid vectors, rmin, 
            and packing fraction.
    """
    from phenum.vector_utils import _minkowski_reduce_basis

    Bs, Hs = find_srBs(A,kpd,exact=True)
    #now we want to find the grid vectors
    B_grids = []
    for B in Bs:
        B_grids.append(np.transpose(np.linalg.inv(B)))
    
    #now I want to find the packing fraction and the r_min for each of these grids.

    grids = []
    for grid in B_grids:
        min_grid = np.transpose(_minkowski_reduce_basis(np.transpose(grid),1E-10))
        rm = max(np.linalg.norm(min_grid,axis=0))
        pf = 4/3.*np.pi*rm**2/np.dot(min_grid[0],np.cross(min_grid[1],min_grid[2]))
        grid={"grid_vecs":grid,"r_min":rm,"packing_frac":pf}
        grids.append(grid)

    return grids
