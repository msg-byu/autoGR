"""Finds the possible symmetry preserving grids for a given lattice.
"""

import numpy as np
from universal import find_srBs
from lat_type import lat_type

def spGrids(A,kpd):
    """Finds the symmetry preserving grids for a given lattice at a given density.

    Args:
        lattice (list of list): The parent lattice vectors as columns of the a matrix.
        kpd (int): The target k-point density.
    
    Returns:
        grids (list of dict): A list of dictionaries contaning the grid vectors, rmin, 
            and packing fraction.
    """
    from phenum.vector_utils import _minkowski_reduce_basis

    Bs = find_srBs(A,kpd,exact=True)
    #now we want to find the grid vectors
    B_grids = []
    ns = []
    for B in Bs:
        B_grids.append(np.transpose(np.linalg.inv(B)))
        ns.append(np.linalg.det(B)/np.linalg.det(A))
    
    #now I want to find the packing fraction and the r_min for each of these grids.

    grids = []
    for grid in B_grids: 
        min_grid = np.transpose(_minkowski_reduce_basis(np.transpose(grid),1E-7))
        rm = min(np.linalg.norm(min_grid,axis=0))
        pf = 4/3.*np.pi*rm**2/np.dot(min_grid[0],np.cross(min_grid[1],min_grid[2]))
        grid={"grid_vecs":min_grid,"r_min":rm,"packing_frac":pf}
        grids.append(grid)

    return grids
