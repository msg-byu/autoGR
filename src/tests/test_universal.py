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
    

def test_find_volumes():
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

    count = 1
    folder = "/Users/wileymorgan/codes/opf_kgrids/src/fortran/tests/find_kgrids/{0}.{1}"

    #hex_12
    A = [[ 0.       ,  1.5      ,  1.       ],
       [ 0.       , -0.8660254,  0.       ],
       [ 2.       ,  2.       ,  2.       ]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_12_14.txt").reshape(14,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #hex_22
    A = [[-0.5      ,  0.5      , -0.5      ],
       [ 0.8660254,  0.8660254,  0.8660254],
       [ 0.       , -0.5      , -0.5      ]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_22_1.txt").reshape(1,3,3)

    assert np.allclose(Bs,true_Bs)
    
    #sc_3
    A = [[0, 1, 1],
       [0, 1, 0],
       [1, 1, 1]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_3_3.txt").reshape(3,3,3)

    assert np.allclose(Bs,true_Bs)
    
    #fcc_1
    A = [[1, 2, 1],
         [1, 2, 2],
         [0, 2, 1]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_1_3.txt").reshape(3,3,3)

    assert np.allclose(Bs,true_Bs)
    
    #bcc_5
    A = [[ 1,  1,  0],
         [ 1,  1,  2],
         [-1,  1,  0]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_5_3.txt").reshape(3,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #rhom_9
    A = [[4, 7, 5],
         [3, 6, 5],
         [3, 7, 5]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_9_1.txt").reshape(1,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #rhom_4_2
    A = [[-1.11652 , -2.11652 , -2.11652 ],
         [-0.610985, -1.933865, -0.610985],
         [ 0.616515, -0.883485, -0.383485]]
    Bs = find_supercells(A,20)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_4_2_13.txt").reshape(13,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #rhom_24
    A = [[-0.255922,  0.255918, -1.255922],
         [-1.44338 , -1.44338 , -1.44338 ],
         [ 0.92259 , -0.922588, -0.07741 ]]
    Bs = find_supercells(A,10)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_24_9.txt").reshape(9,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #stet_11
    A = [[0, 1, 1],
         [0, 1, 0],
         [2, 2, 2]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_11_3.txt").reshape(3,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #stet_21
    A = [[ 0. ,  1. ,  0. ],
         [ 1. ,  1. ,  1. ],
         [ 0. ,  0.5,  0.5]]
    Bs = find_supercells(A,20)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_21_20.txt").reshape(20,3,3)

    assert np.allclose(Bs,true_Bs)
        
    #body_tet_15
    A = [[ 1,  1,  0],
         [ 1,  1,  2],
         [-2,  2,  0]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_15_1.txt").reshape(1,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #body_tet_7
    A = [[ 1.95095 ,  1.      ,  0.      ],
         [ 1.19163 ,  1.60788 ,  2.60788 ],
         [ 0.879663, -1.55394 ,  0.44606 ]]
    Bs = find_supercells(A,20)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_7_15.txt").reshape(15,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #body_tet_6
    A = [[ 1.80278 ,  1.80278 ,  0.80278 ],
         [-1.47253 ,  1.13535 , -0.47253 ],
         [ 0.762655,  1.208715,  2.762655]]
    Bs = find_supercells(A,20)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_6_15.txt").reshape(15,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #body_tet_18
    A = [[-2, -1, -2],
         [-1, -3, -1],
         [ 1,  4,  3]]
    Bs = find_supercells(A,20)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_18_15.txt").reshape(15,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #so_32
    A = [[0, 1, 1],
         [0, 2, 0],
         [3, 3, 3]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_32_21.txt").reshape(21,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #so_32
    A = [[0, 1, 1],
         [0, 2, 0],
         [3, 3, 3]]
    Bs = find_supercells(A,5)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_32_70.txt").reshape(70,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #base_ortho_38_13
    A = [[ 0. ,  1. ,  0.5],
         [ 0. ,  0. ,  1. ],
         [ 3. ,  3. ,  3. ]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_38_13_9.txt").reshape(9,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #base_ortho_23
    A = [[ 2.       ,  2.6666667,  1.6666667],
         [-1.       , -1.54116  , -2.54116  ],
         [-1.       ,  1.87449  ,  0.87449  ]]
    Bs = find_supercells(A,5)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_23_42.txt").reshape(42,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #base_ortho_40
    A = [[-1.05557 ,  1.56246 , -0.05557 ],
         [ 1.99895 ,  2.380916,  2.99895 ],
         [-0.943376, -0.943376,  0.056624]]
    Bs = find_supercells(A,5)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_40_42.txt").reshape(42,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #base_ortho_36
    A = [[-1.43541,  0.9788 , -0.43541],
         [-1.43541, -1.84962, -0.43541],
         [ 1.37083,  2.37083,  2.37083]]
    Bs = find_supercells(A,5)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_36_42.txt").reshape(42,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #body_ortho_19
    A = [[ 0. ,  0.5,  0.5],
         [ 0. ,  3. ,  1. ],
         [ 3. ,  4.5,  4.5]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_19_3.txt").reshape(3,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #body_ortho_8
    A = [[ 3.41012  ,  3.82288  ,  4.82156  ],
         [-1.1237578,  1.1771244, -1.0351956],
         [-1.268178 , -2.       , -3.268178 ]]
    Bs = find_supercells(A,5)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_8_34.txt").reshape(34,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #body_ortho_42
    A = [[ 1.61803,  1.0817 ,  0.0817 ],
         [-0.61803,  1.74903,  0.74903],
         [-1.     , -1.33073, -2.33073]]
    Bs = find_supercells(A,5)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_42_34.txt").reshape(34,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #face_ortho_26
    A = [[ 0. ,  0.5,  0. ],
         [ 0. ,  1. ,  1. ],
         [ 3. ,  6. ,  4.5]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_26_9.txt").reshape(9,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #face_ortho_16
    A = [[ 1.04442 ,  1.515172,  2.824216],
         [ 1.43973 , -2.79933 ,  1.25993 ],
         [ 1.68415 , -0.68415 ,  1.68415 ]]
    Bs = find_supercells(A,5)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_16_34.txt").reshape(34,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #simple_mono 33
    A = [[ 2. , -1.5,  2.5],
         [ 0. ,  2. ,  0. ],
         [ 0. ,  2. ,  2. ]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_33_49.txt").reshape(49,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #simple_mono 35
    A = [[ 1.      , -0.050882,  0.331088],
         [ 1.      ,  0.348726,  2.96676 ],
         [ 1.      , -3.29785 , -0.29785 ]]
    Bs = find_supercells(A,10,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_35_49.txt").reshape(49,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #simple_mono 34
    A = [[ 1.        ,  0.05875978,  0.83401491],
         [ 1.        , -3.86782784, -0.64308297],
         [ 1.        , -0.19093194,  2.80906806]]
    Bs = find_supercells(A,1)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_34_77.txt").reshape(77,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #base_mono 14
    A = [[ 1. , -0.5,  1.5],
         [ 1. ,  1. ,  1. ],
         [ 0. ,  2. ,  2. ]]
    Bs = find_supercells(A,8,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_14_31.txt").reshape(31,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #base_mono 10 17
    A = [[-1.46391,  2.46391, -1.46391],
         [ 0.     ,  3.     ,  2.     ],
         [ 1.96391, -0.96391,  1.96391]]
    Bs = find_supercells(A,8,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_10_17_31.txt").reshape(31,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #base_mono 20 25
    A = [[ 1.       ,  1.39898  ,  1.69779  ],
         [ 1.       , -3.8834405, -0.4322505],
         [ 1.       ,  3.23446  ,  4.23446  ]]
    Bs = find_supercells(A,8,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_20_25_31.txt").reshape(31,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #base_mono 27
    A = [[ 0.464824, -3.082857, -0.535176],
         [-1.464824,  1.082857, -2.464824],
         [-1.907413,  0.907413, -1.907413]]
    Bs = find_supercells(A,8,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_27_31.txt").reshape(31,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #base_mono 28
    A = [[-1.44896 ,  0.791384, -1.106536],
         [ 0.948958, -3.291382, -0.393466],
         [-1.      , -1.02006 , -3.02006 ]]
    Bs = find_supercells(A,8,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_28_31.txt").reshape(31,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #base_mono 29 30
    A = [[-0.666125,  3.284155,  0.951905],
         [ 1.16613 , -0.784164,  0.548096],
         [ 2.04852 , -1.04852 ,  3.04852 ]]
    Bs = find_supercells(A,8,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_29_30_31.txt").reshape(31,3,3)

    assert np.allclose(Bs,true_Bs)
                               
    #base_mono 37 39 and 41
    A = [[-1.79092,   1.,        1.      ],
        [-1.47209,  -1.41421,   0.      ],
        [ 0.790922, -1.,        1.      ]]
    Bs = find_supercells(A,8,exact=True)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_37_39_41_31.txt").reshape(31,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #base_mono 43
    A = [[-0.39716,  1.64194, -1.79432],
         [-0.34718, -1.14194, -1.69436],
         [ 2.49434, -1.     ,  3.98868]]
    Bs = find_supercells(A,1)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_43_41.txt").reshape(41,3,3)

    assert np.allclose(Bs,true_Bs)
                
    #triclinic
    A = [[ 2.24478,  3.48956, -1.79432],
         [-0.48912, -1.97824, -1.69436],
         [ 2.49434,  3.98868,  3.98868]]
    Bs = find_supercells(A,5)

    true_Bs = np.loadtxt("tests/test_output/find_supercells_basis_44_31.txt").reshape(31,3,3)

    assert np.allclose(Bs,true_Bs)
