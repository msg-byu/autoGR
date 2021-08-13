import numpy as np
import pytest

from opf_python.pyniggli import reduced_cell

def test_reduction():
    """Tests the niggli reduction of some basic cells."""

    i = 1
    import csv
    
    #1
    A = np.transpose([[1,0,0],[0,1,0],[0,0,2]])
    B = reduced_cell(A)
    assert np.allclose(B.C,np.array([[1,0,0],[0,1,0],[0,0,1]]))
    assert np.allclose(B.niggli,A)
    assert np.allclose(B.niggli,np.dot(A,B.C))
    
    #2
    A = np.transpose([[0.5,0,0.5],[0,3,0],[0.5,0,-0.5]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli,np.transpose([[-0.5,0,-0.5],[-0.5,0,0.5],[0,-3,0]]))
    assert np.allclose(B.niggli,np.dot(A,B.C))

    #3
    A = np.transpose([[1.00000000,0.00000000,0.00000000],[-0.50000000,0.86602540,1.63299320],[0.00000000,-1.73205080,1.63299320]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli,[[-1.       ,  0.5      ,  0.       ],
       [ 0.       , -0.8660254, -1.7320508],
       [ 0.       , -1.6329932,  1.6329932]])

    #4
    A = np.transpose([[1.00000000, 0.00000000, 0.00000000],
                      [-0.50000000, 0.86602540, 0.00000000],
                      [0.00000000, 0.00000000, 3.26598640]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli, np.transpose([[ 1., 0., 0.],
                                               [-0.5, 0.8660254, 0.],
                                               [ 0., 0., 3.2659864]]))

    #5
    A = np.transpose([[1.00000000, 0.00000000, 0.00000000],
                      [-0.50000000, 0.86602540, 1.63299320],
                      [0.00000000, -1.73205080, 1.63299320]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli, np.transpose([[-1., 0., 0.],
                                               [0.5, -0.8660254, -1.6329932],
                                               [0., -1.7320508, 1.6329932]]))

    #6
    A = np.transpose([[1.00000000, 0.00000000, 0.00000000],
                      [0.50000000, -0.86602540, 3.26598640],
                      [0.00000000, -1.73205080, 0.00000000]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli, np.transpose([[1., 0., 0.],
                                               [0., -1.7320508, 0.],
                                               [-0.5, 0.8660254, -3.2659864]]))

    #7
    A = np.transpose([[1.00000000, 0.00000000, 0.00000000],
                      [0.50000000, 4.33012700, 0.00000000],
                      [0.00000000, 0.00000000, 1.63299320]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli, np.transpose([[-1., 0., 0.],
                                               [0., 0., 1.6329932],
                                               [0.5, 4.330127, 0.]]))

    #8
    A = np.transpose([[1.00000000, 0.00000000, 0.00000000],
                      [0.50000000, -0.86602540, 4.89897960],
                      [0.00000000, -1.73205080, 0.00000000]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli, np.transpose([[1., 0., 0.],
                                               [0., -1.7320508, 0.],
                                               [-0.5, 0.8660254, -4.8989796]]))

    #9
    A = np.transpose([[0.00000000, -1.73205080, 1.63299320],
                      [0.50000000, 2.59807620, 3.26598640],
                      [1.00000000, 0.00000000, 0.00000000]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli, np.transpose([[-1., 0., 0.],
                                               [0., 1.7320508, -1.6329932],
                                               [0.5, 2.5980762, 3.2659864]]))

    #10
    A = np.transpose([[0.5, 0.5, -0.5],
                      [-0.5, 0.5, 0.5],
                      [1.0, 0.0, 1.0]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli, np.transpose([[0.5, 0.5, -0.5],
                                               [-0.5, 0.5, 0.5],
                                               [1., 0., 1.]]))

    #11
    A = np.transpose([[1.00000000, 0.00000000, 0.00000000],
                      [0.00000000, 0.00000000, 1.00000000],
                      [0.50000000, -1.50000000, 0.50000000]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli, np.transpose([[-1., 0., 0.],
                                               [0., 0., -1.],
                                               [0.5, -1.5, 0.5]]))

    #12
    A = np.transpose([[.05, 2.7, 3.3],
                      [0.1, 0.7, 4.5],
                      [.99, .3, 5.4]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli, np.transpose([[-0.89,  0.4 , -0.9 ],
                                               [-0.84, -1.6 ,  0.3 ],
                                               [-1.68,  1.5 ,  2.7 ]]))

    #13
    A = np.transpose([[1, -.1, 0],
                      [-0.3, 1, .3],
                      [-.3, -0.1, -1.5 ]])
    B = reduced_cell(A)
    assert np.allclose(B.niggli, np.transpose([[-1., 0.1, 0.],
                                               [0.3, -1., -0.3],
                                               [0.4, 0.8, -1.2]]))

    #14
    A = np.array([[ 0.03682244,  0.        ,  0.        ],
       [ 0.        ,  0.06377834,  0.        ],
       [ 0.96209276,  0.96209276,  1.92418552]])
    B = reduced_cell(A,eps=1E-10)
    assert np.allclose(B.niggli, [[-0.07364488, -0.03682244, -0.03682244],
                                  [ 0.        ,  0.06377834,  0.        ],
                                  [ 0.        ,  0.        , -0.96209276]])

    # 15
    A = np.array([[ 0.07364488,  0.        ,  0.        ],
                  [ 0.        ,  0.12755668,  0.        ],
                  [ 0.24052319,  0.24052319,  0.48104638]])
    B = reduced_cell(A,eps=1E-10)
    assert np.allclose(B.niggli, [[-0.14728977, -0.07364488, -0.07364488],
                                  [ 0.        ,  0.12755668,  0.        ],
                                  [ 0.        ,  0.        , -0.24052319]])

    #16
    A = np.array([[ 1.20559446,  0.        ,  0.        ],
                  [ 0.        ,  2.08815085,  0.        ],
                  [ 3.93745511,  3.93745511,  7.87491022]])
    B = reduced_cell(A,eps=1E-5)
    assert np.allclose(B.niggli,[[-2.41118892, -1.20559446, -1.20559446],
                                 [ 0.        ,  2.08815085,  0.        ],
                                 [ 0.        ,  0.        , -3.93745511]])
    
    with pytest.raises(ValueError):
        reduced_cell([[0,0,0],[0,0,0],[0,0,0]])
    with pytest.raises(ValueError):
        reduced_cell([0,0,0])
    with pytest.raises(ValueError):
        reduced_cell([[2,2,2],[1,1,1]])

def test_findC3():
    
    A = np.array([[1,0,0],[0,1,0],[0,0,1]])

    B = reduced_cell(A,eps=1E-7)

    C = B._find_C3(-1,-1,-1)
    assert np.allclose(C,np.array([[-1,0,0],[0,-1,0],[0,0,-1]]))

    C = B._find_C3(1,1,1)
    assert np.allclose(C,np.array([[1,0,0],[0,1,0],[0,0,1]]))

def test_findC4():
    
    A = np.array([[1,0,0],[0,1,0],[0,0,1]])

    B = reduced_cell(A)

    C = B._find_C4(-1,-1,-1)
    assert np.allclose(C,np.array([[1,0,0],[0,1,0],[0,0,1]]))
    
    C = B._find_C4(0,0,1)
    assert np.allclose(C,np.array([[1,0,0],[0,-1,0],[0,0,-1]]))
    
    C = B._find_C4(1,-1,0)
    assert np.allclose(C,np.array([[-1,0,0],[0,1,0],[0,0,-1]]))
    
    C = B._find_C4(0,1,-1)
    assert np.allclose(C,np.array([[-1,0,0],[0,-1,0],[0,0,1]]))

def test_niggli_check():

    A = np.array([[1,0,0],[0,1,0],[0,0,1]])

    Bc = reduced_cell(A,path_=True)

    eps = 1E-3

    A, B, C, xi, eta, zeta = 0, 2, 1, 0, 0, 0
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 1, 2, 1, 0, 0
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 2, 2, 0, 2, 0
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 1, 2, 0, 1, 0
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 2, 3, -3, -1, 0
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 2, 3, -1, -3, -2
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 2, 2.1, -1.9, -.9, -.9
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 2, 2.1, 2, .1, .9
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 2, 2.1, .1, 1, .9
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 2, 2.1, .1, .5, 1
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 2, 2.1, -2, -.1, -.9
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 2, 2.1, -.1, -1, -.9
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 2, 2.1, -.1, -.5, -1
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
    A, B, C, xi, eta, zeta = 1, 2, 2.1, -1.9, -.8, -.3
    assert not Bc._niggli_check(A,B,C,xi,eta,zeta,eps)
    
