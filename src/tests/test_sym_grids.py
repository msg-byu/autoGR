"""Tests the srHNF code for each Bravais Lattice."""

import pytest
import sys
import numpy as np

def test_btet_grid():
    """Tests the grid generation of the grids."""
    from opf_python.sym_grids import spGrids
    
    a1 = np.array([-1,1,2])
    a2 = np.array([1,-1,2])
    a3 = np.array([1,1,-2])

    lat_param = np.random.rand()*10
    A = lat_param*np.transpose([a1,a2+a3,a1+a3])
    B = np.linalg.inv(np.transpose(A))
    kpd = np.random.randint(10,high=100)
    grids = spGrids(A,kpd)

    for grid in grids:
        vecs = grid["grid_vecs"]
        test = np.dot(np.linalg.inv(vecs),B)
        assert np.allclose(abs(test-np.round(test)),0)
        assert np.allclose(abs(np.linalg.det(B)/np.linalg.det(vecs)),kpd)

def test_fcc_grid():
    """Tests the grid generation of the grids."""
    from opf_python.sym_grids import spGrids
    
    lat_param = np.random.rand()*10
    A = lat_param*np.transpose([[2,2,2],[1,1,2],[1,2,1]])
    B = np.transpose(np.linalg.inv(A))

    kpd = np.random.choice([4,16,8,32,128,27,108,432])
    grids = []

    grids = spGrids(A,kpd)

    for grid in grids:
        vecs = grid["grid_vecs"]
        test = np.dot(np.linalg.inv(vecs),B)
        assert np.allclose(abs(test-np.round(test)),0)
        assert np.allclose(abs(np.linalg.det(B)/np.linalg.det(vecs)),kpd)
        
def test_bcc_grid():
    """Tests the grid generation of the grids."""
    from opf_python.sym_grids import spGrids
    
    lat_param = np.random.rand()*10
    A = lat_param*np.transpose([[1,-1,3],[1,-1,1],[1,1,-1]])
    B = np.linalg.inv(np.transpose(A))

    kpd = np.random.choice([2,4,16,8,32,54,27,108])

    grids = spGrids(A,kpd)

    for grid in grids:
        vecs = grid["grid_vecs"]
        test = np.dot(np.linalg.inv(vecs),B)
        assert np.allclose(abs(test-np.round(test)),0)        
        assert np.allclose(abs(np.linalg.det(B)/np.linalg.det(vecs)),kpd)
