import numpy as np
import pytest

def test_transform_HNFs():
    """Tests the transformation of the symmetry preserving HNFs from our
    basis to the niggli basis."""
    from opf_python.universal import transform_supercells

    # sc case (no actual transformation).
    Nu = np.array([[1,0,0],[0,1,0],[0,0,1]])
    Cu= np.array([[1,0,0],[0,1,0],[0,0,1]])
    Co = np.array([[1,0,0],[0,1,0],[0,0,1]])
    No = np.array([[1,0,0],[0,1,0],[0,0,1]])
    O = np.array([[1,0,0],[0,1,0],[0,0,1]])

    spHNFs = np.array([[[1,0,0],[0,1,0],[1,1,2]],[[1,0,0],[1,2,0],[1,0,2]]])

    new = transform_supercells(spHNFs,No,Nu,Co,Cu,O)

    assert np.allclose(new,spHNFs)
    
    #sc case (with transformation).
    Nu = np.array([[-1,0,0],[0,0,1],[0,1,0]])
    Cu= np.array([[-1,0,-1],[0,0,1],[0,1,0]])
    No = np.array([[1,0,0],[0,1,0],[0,0,1]])
    Co = np.array([[1,0,0],[0,1,0],[0,0,1]])
    O = np.array([[1,0,0],[0,1,0],[0,0,1]])

    spHNFs = np.array([[[1,0,0],[0,1,0],[1,1,2]],[[1,0,0],[1,2,0],[1,0,2]]])

    new = transform_supercells(spHNFs,No,Nu,Co,Cu,O)

    bus = np.array([[[1,1,0],[-1,1,1],[0,0,1]],[[1,1,0],[-1,1,0],[-1,-1,2]]])

    assert np.allclose(new,bus)

    # monoclinic case.
    Nu = np.array([[-1,-1,0.5],[-1,1,0],[0,0,2]])
    Cu= np.array([[0,1,0],[-3,-9,-4],[2,6,3]])
    No = np.array([[1,1,-0.5],[1,-1,0],[0,0,-2]])
    Co = np.array([[1,1,0],[0,-1,0],[0,0,-1]])
    O = np.array([[1,0,0.5],[1,2,0],[0,0,2]])

    spHNFs = np.array([[[1,0,0],[0,1,0],[1,0,4]],
                       [[1,0,0],[0,1,0],[2,0,4]],
                       [[1,0,0],[0,1,0],[3,0,4]],
                       [[2,0,0],[0,1,0],[0,0,2]],
                       [[2,0,0],[0,1,0],[0,1,2]],
                       [[2,0,0],[0,1,0],[1,0,2]],
                       [[2,0,0],[0,1,0],[1,1,2]],
                       [[2,0,0],[0,2,0],[0,0,1]],
                       [[2,0,0],[1,2,0],[0,0,1]]])

    new = transform_supercells(spHNFs,No,Nu,Co,Cu,O)

    bus = np.array([[[3,8.5,12],[4,3,4],[4,22,32]],
                    [[4,10,14],[4,3,4],[8,28,40]],
                    [[5,11.5,16],[4,3,4],[12,34,48]],
                    [[4,8,11],[6,6,8],[0,8,12]],
                    [[4.5,8,11],[6,6,8],[2,8,12]],
                    [[5,9.5,13],[6,6,8],[4,14,20]],
                    [[5.5,9.5,13],[6,6,8],[6,14,20]],
                    [[4,7,9.5],[8,6,8],[0,4,6]],
                    [[4,7,9.5],[12,12,16],[0,4,6]]])

    assert np.allclose(new,bus)
    

def test_find_volemus():
    """Tests the transformation of the symmetry preserving HNFs from our
    basis to the niggli basis."""
    from opf_python.universal import find_volumes

    a,b = find_volumes(1,10)
    assert np.allclose(a,[1, 4, 16, 8, 32, 128, 27, 108, 432, 64, 256, 125, 500, 216, 864, 343, 512, 729, 1000])
    assert b is None 

    a,b = find_volumes(3,10)
    assert np.allclose(a,[1, 2, 4, 8, 16, 32, 27, 54, 108, 64, 128, 256, 125, 250, 500, 216, 432, 864, 343, 686, 512, 729, 1000])
    assert b is None 


    a,b = find_volumes(31,111)
    assert a==[13]
    assert b==2

    a,b = find_volumes(31,1110)
    assert a==[41]
    assert b==3

    a,b = find_volumes(10,50)
    assert np.allclose(a,[45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55])
    assert b is None

    a,b = find_volumes(20,4)
    assert np.allclose(a,[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13])
    assert b is None

def test_find_supers():
    """Finds the supercells of the given cell."""

    from opf_python.universal import find_supercells

    #hex_12
    A = [[ 0.       ,  1.5      ,  1.       ],
       [ 0.       , -0.8660254,  0.       ],
       [ 2.       ,  2.       ,  2.       ]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_hex_12_13_3_3.txt").reshape(13,3,3)

    assert np.allclose(Bs,true_Bs)

    #hex_22
    A = [[-0.5      ,  0.5      , -0.5      ],
       [ 0.8660254,  0.8660254,  0.8660254],
       [ 0.       , -0.5      , -0.5      ]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_hex_22_1_3_3.txt").reshape(1,3,3)

    assert np.allclose(Bs,true_Bs)
    
    #sc_3
    A = [[0, 1, 1],
       [0, 1, 0],
       [1, 1, 1]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_sc_3_23.txt").reshape(23,3,3)

    assert np.allclose(Bs,true_Bs)
    
