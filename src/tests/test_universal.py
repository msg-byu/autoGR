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
    assert np.allclose(a,[16, 27, 32])
    assert b is None 

    a,b = find_volumes(3,10)
    assert np.allclose(a,[16, 27, 32])
    assert b is None 


    a,b = find_volumes(31,111)
    assert a==[13]
    assert b==2

    a,b = find_volumes(31,1110)
    assert a==[41]
    assert b==3

    a,b = find_volumes(10,50)
    assert np.allclose(a,[50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60])
    assert b is None

def test_find_supers():
    """Finds the supercells of the given cell."""

    from opf_python.universal import find_supercells

    #hex_12
    A = [[ 0.       ,  1.5      ,  1.       ],
       [ 0.       , -0.8660254,  0.       ],
       [ 2.       ,  2.       ,  2.       ]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_hex_12_14.txt").reshape(14,3,3)

    assert np.allclose(Bs,true_Bs)

    #hex_22
    A = [[-0.5      ,  0.5      , -0.5      ],
       [ 0.8660254,  0.8660254,  0.8660254],
       [ 0.       , -0.5      , -0.5      ]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_hex_22_1.txt").reshape(1,3,3)

    assert np.allclose(Bs,true_Bs)
    
    #sc_3
    A = [[0, 1, 1],
       [0, 1, 0],
       [1, 1, 1]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_sc_3_3.txt").reshape(3,3,3)

    assert np.allclose(Bs,true_Bs)
    
    #fcc_1
    A = [[1, 2, 1],
         [1, 2, 2],
         [0, 2, 1]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_fcc_1_3.txt").reshape(3,3,3)

    assert np.allclose(Bs,true_Bs)
    
    #bcc_5
    A = [[ 1,  1,  0],
         [ 1,  1,  2],
         [-1,  1,  0]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_bcc_5_3.txt").reshape(3,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #rhom_9
    A = [[4, 7, 5],
         [3, 6, 5],
         [3, 7, 5]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_rhom_9_1.txt").reshape(1,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #rhom_4_2
    A = [[-1.11652 , -2.11652 , -2.11652 ],
         [-0.610985, -1.933865, -0.610985],
         [ 0.616515, -0.883485, -0.383485]]
    Bs = find_supercells(A,20)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_rhom_4_2_13.txt").reshape(13,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #rhom_24
    A = [[-0.255922,  0.255918, -1.255922],
         [-1.44338 , -1.44338 , -1.44338 ],
         [ 0.92259 , -0.922588, -0.07741 ]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_rhom_24_9.txt").reshape(9,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #stet_11
    A = [[0, 1, 1],
         [0, 1, 0],
         [2, 2, 2]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_stet_11_3.txt").reshape(3,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #stet_21
    A = [[ 0. ,  1. ,  0. ],
         [ 1. ,  1. ,  1. ],
         [ 0. ,  0.5,  0.5]]
    Bs = find_supercells(A,20)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_stet_21_20.txt").reshape(20,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #body_tet_15
    A = [[ 1,  1,  0],
         [ 1,  1,  2],
         [-2,  2,  0]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_body_tet_15_1.txt").reshape(1,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #body_tet_7
    A = [[ 1.95095 ,  1.      ,  0.      ],
         [ 1.19163 ,  1.60788 ,  2.60788 ],
         [ 0.879663, -1.55394 ,  0.44606 ]]
    Bs = find_supercells(A,20)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_body_tet_7_15.txt").reshape(15,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #body_tet_6
    A = [[ 1.80278 ,  1.80278 ,  0.80278 ],
         [-1.47253 ,  1.13535 , -0.47253 ],
         [ 0.762655,  1.208715,  2.762655]]
    Bs = find_supercells(A,20)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_body_tet_6_15.txt").reshape(15,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #body_tet_18
    A = [[-2, -1, -2],
         [-1, -3, -1],
         [ 1,  4,  3]]
    Bs = find_supercells(A,20)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_body_tet_18_15.txt").reshape(15,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #so_32
    A = [[0, 1, 1],
         [0, 2, 0],
         [3, 3, 3]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_so_32_21.txt").reshape(21,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #so_32
    A = [[0, 1, 1],
         [0, 2, 0],
         [3, 3, 3]]
    Bs = find_supercells(A,5)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_so_32_70.txt").reshape(70,3,3)

    assert np.allclose(Bs,true_Bs)
                
