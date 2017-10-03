import numpy as np
import pytest
from opf_python.niggli_lat_id import niggli_id

def test_bo():
    """Tests of body centered orthorhombic lattice identiication."""

    A = np.transpose([[-0.5,1,1.5],[0.5,-1,1.5],[0.5,1,-1.5]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 19

    A = np.transpose([[1, 1, 2],[1.41144,0.0885622,-2],
                 [-1.99868,1.21232,-0.731822]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 8

    A = np.transpose([[1, 1, 1],[1.61803,-0.61803,-1],
                 [-1.53633,1.36706,-1.33073]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 42

    # Test of failed case.
    A = [[ 0.07364488,  0.        ,  0.        ],
         [ 0.        ,  0.06377834,  0.        ],
       [ 0.48104638,  0.48104638,  0.96209276]]
    lat_t, nig_n, lat_f = niggli_id(A,eps_=1E-10)
    assert lat_t == 'face centered orthorhombic'
    assert lat_f == 5

def test_fcc():
    """Tests of fcc lattice identification."""

    A = np.transpose([[0.5,0,0.5],[0.5,0.5,0],[0,0,1]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 1

def test_sc():
    """Tests of sc lattice identification."""

    A = np.transpose([[1,0,0],[1,1,0],[0,0,1]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 3
    
def test_bcc():
    """Tests of bcc lattice identification."""

    A = [[-1,1,1],[1,-1,1],[1,1,-1]]
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 5

def test_stet():
    """Tests of simple tet lattice identification."""

    a1 = np.array([1,0,0])
    a2 = np.array([0,1,0])
    a3 = np.array([0,0,2])

    As = [[a1+a2,a2,a3+a1],[a1+a2+a3,a2+a3,a1+a2],[a1,a2+a1,a3+a3],[-a1,-a2+a3,a3]]
    As = [np.transpose(i) for i in As]

    for A in As:
        lat_t, nig_n, lat_f = niggli_id(A)
        assert nig_n == 11

    A = np.transpose([[0,0,0.5],[1,0,0],[0,1,0]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 21

def test_so():
    """Tests of simple orthorhombic lattice identification."""

    As = [[[1.00000000, 0.00000000, 0.00000000],
       [0.00000000, 1.73205080, 0.00000000],
       [0.00000000, 0.00000000, 3.26598640]],
      [[0.50000000, 0.50000000, 0.00000000],
      [0.00000000, 0.00000000, 1.00000000],
      [1.50000000, -1.50000000, 0.00000000]],
      [[0.50000000, 0.50000000, 0.00000000],
      [0.00000000, 0.00000000, 1.00000000],
      [2.00000000, -2.00000000, 0.00000000]]]
    As = [np.transpose(i) for i in As]
    
    for A in As:
        lat_t, nig_n, lat_f = niggli_id(A)
        assert nig_n == 32

def test_rhom():
    """Tests of rhombohedral lattice identification."""

    A = np.transpose([[1,2,2],[2,1,2],[4,3,3]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 9

    A = np.transpose([[-1,0,-1],[1.51184,0,-0.845178],
                      [-0.255922,-1.44338,0.92259]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 24
     
    A = np.transpose([[-1, 0,-1],[0, -1.32288, -0.5],
                      [-1.11652, -0.610985, 0.616515]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 2
    
    A = np.transpose([[-1,0,-1],[0, -1.32288, 0.5],
                      [-0.548584, 0.774292, 1.04858]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 4

def test_hex():

    A = np.transpose([[1,0,0],[0.5,-0.8660254037844386,0],
                      [0,0,2]])
    lat_t, nig_n, lat_f = niggli_id(A,G=[1.0, 0.99999999999999989, 4.0, 0.0, 0.0, -0.5])
    assert nig_n == 12

    A = np.transpose([[0,0,-0.5],[1,0,0],
                      [-0.5,0.8660254037844386,0]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 22
    
def test_bctet():
    """Tests of body centered tetragonal lattice identification."""
    A = np.transpose([[-1,1,2],[1,-1,2],[1,1,-2]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 15

    A = np.transpose([[-1, 1, 2],[1, 1.60788, -1.55394],
                     [1.95095, -1.41625, 0.433603]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 7

    A = np.transpose([[-1, 1, 2],[1, 1.60788, -1.55394],
                     [1.80278, -1.47253, 0.762655]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 6

    A = np.transpose([[0,0,2],[1,-2,1],
                  [-2,-1,1]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 18

def test_fco():
    """Tests of face centered orthorhombic lattice identification."""

    A = np.transpose([[0,1,1.5],[0.5,0,1.5],[0,0,3]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 26

    A = np.transpose([[1,1,-1],[0.779796,-1.1798,1],
                 [-1.04442,-1.43973,-1.68415]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 16

def test_cco():
    """Tests of base centered orthorhombic lattice identification."""

    A = np.transpose([[0.5,1,0],[0.5,-1,0],[0,0,3]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 38

    A = np.transpose([[1,1,1],[1,-1,-1],
                      [0,-1.73205,1.73205]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 13

    A = np.transpose([[1,1,1],[2,-1,-1],
                      [-0.3333333,-1.54116,1.87449]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 23

    A = np.transpose([[1,1,1],
                      [1.61803,-0.618034,-1],
                      [-1.05557,1.99895,-0.943376]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 40

    A = np.transpose([[1, 1, 1],
                      [1.41421, -1.41421, 0],
                      [-1.43541, -1.43541, 1.37083]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 36
    
def test_sm():
    """Tests of the simple monoclinic lattice identification."""

    A = np.transpose([[2,0,0],
                      [0,2,0],
                      [0.5,0,2]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 33

    A = np.transpose([[1,1,1],
                      [1.61803,-0.618034,-1],
                      [-0.668912,1.96676,-1.29785]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 35

    A = np.transpose([[1,1,1],
                      [1.22474487,-1.22474487,-1],
                      [-0.165985087,-1.64308297,1.80906806]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 34

def test_cm():
    """Tests of the base centered monoclinic lattice identification."""

    A = np.transpose([[1,1,0],
                      [1,-1,0],
                      [0.5,0,2]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 14

    A = np.transpose([[1,1,1],
                  [1,-1,1],
                  [-1.46391,0,1.96391]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 10

    A = np.transpose([[1,1,1],
                      [0.809568,-1.15957,-1],
                      [-1.05387,-1.61088,1.51474]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 17

    A = np.transpose([[1,1,1],
                  [1.70119,-1.45119,1],
                  [-1.0034,0.0189395,2.23446]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 20

    A = np.transpose([[1,1,1],
                      [1.45119,-1.70119,-1],
                      [-1.16241,-1.56776,1.48018]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 25

    A = np.transpose([[1,1,0],
                      [1.618033,-0.618033,1],
                      [-0.464824,1.464824,1.907413]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 27

    A = np.transpose([[1,1,0],
                      [1.44896,-0.948958,1],
                      [-0.342424,1.342424,2.02006]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 28

    A = np.transpose([[1,1,0],
                      [1.61803,-0.618034,1],
                      [-0.666125,1.16613,2.04852]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 29

    A = np.transpose([[1,1,0],
                      [1.61803,-0.618034,1],
                      [-0.0361373,0.536137,2.38982]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 30

    A = np.transpose([[1,0,1],
                      [1,-1.41421,-1],
                      [-1.85397,-0.854143,1.35397]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 41

    A = np.transpose([[1,0,1],
                      [1,-1.41421,-1],
                      [-1.79092,-1.47209,0.790922]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 37

    A = np.transpose([[1,0,1],
                      [0,-1.73205,-1],
                      [-1.66542,-0.672857,1.66542]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 39

    A = np.transpose([[1,1,1],
                      [1.64194,-1.14194,-1],
                      [-1.39716,-1.34718,1.49434]])
    lat_t, nig_n, lat_f = niggli_id(A)
    assert nig_n == 43

def test_tric():
    """Test the triclinic littice identification."""

    A = [[1,0,0],[0,1,0],[0,0,1]]
    lat_t, nig_n, lat_f = niggli_id(A,G=[1, 2, 3, .1, .2, .3])
    assert nig_n == 31
    
    A = [[1,0,0],[0,1,0],[0,0,1]]
    lat_t, nig_n, lat_f = niggli_id(A,G=[1, 2, 3, -.1, -.2, -.3])
    assert nig_n == 44

