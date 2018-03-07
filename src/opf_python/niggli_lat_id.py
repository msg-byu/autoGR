"""The methods contained in this module are used to identify a lattice
using the niggli reduced cell."""

import numpy as np

def niggli_id(lattice,eps_=None,G=None):
    """Returns the lattice type, niggli cell number, and crystal family
    number, .i.e, cubic -> 1, hexagonal -> 2, rhombehedral -> 3,
    tetragonal -> 4, orthorhombic -> 5, monoclinic -> 6, and triclinic
    -> 7. 

    Args:
        lattice (numpy array): The lattice vectors as columns of a matrix.
        eps_ (float, optional): The floating point tollerance for comparisons.

    Returns:
        lat_type, niggli_num, lat_fam, basis (str, int, int, list): The lattice type, the 
            niggli case number from the 'International Tables of Crystalography',
            the lattice family type, and the 'canonical' basis choice.
    """

    from opf_python.pyniggli import reduced_cell

    if not isinstance(lattice,np.ndarray):
        lattice = np.array(lattice)

    if eps_ is None:
        eps = 1E-5
    else:
        eps = eps_

    if G is None:
        if eps_ is None:
            temp = reduced_cell(lattice)
        else: 
            temp = reduced_cell(lattice, eps = eps)
        niggli = temp.niggli

        a = niggli[:,0]
        b = niggli[:,1]
        c = niggli[:,2]

        A = np.dot(a,a)
        B = np.dot(b,b)
        C = np.dot(c,c)
        D = np.dot(b,c)
        E = np.dot(a,c)
        F = np.dot(a,b)
    else:
        A = G[0]
        B = G[1]
        C = G[2]
        D = G[3]
        E = G[4]
        F = G[5]

    positive = False
    if D-eps > 0 and E-eps > 0 and F-eps > 0:
        positive = True

    niggli_num = None
    if np.allclose(A,B,atol=eps) and np.allclose(B,C,atol=eps):
        if positive:
            if np.allclose(D,E,atol=eps) and np.allclose(D,F,atol=eps):
                if np.allclose(A/2., D,atol=eps) and np.allclose(A/2., E,atol=eps) and np.allclose(A/2., F,atol=eps):
                    lat_fam = 1
                    lat_type = 'face centered cubic'
                    niggli_num = 1
                    basis = np.transpose([[0,1,1],[1,0,1],[1,1,0]])
                else:
                    lat_fam = 3
                    lat_type = 'rhombohedral'
                    niggli_num = 2
                    basis = np.transpose([[-1, 0, -1],[0, -1.32288, -0.5],
                                          [-1.11652, -0.610985, 0.616515]])
        else:
            if np.allclose(D,E,atol=eps) and np.allclose(D,F,atol=eps):
                if np.allclose(0,D,atol=eps):
                    lat_fam = 1
                    lat_type = 'simple cubic'
                    niggli_num = 3
                    basis = np.transpose([[1,0,0],[0,1,0],[0,0,1]])
                elif np.allclose(-A/3., D,atol=eps):
                    lat_fam = 1
                    lat_type = 'body centered cubic'
                    niggli_num = 5
                    basis = np.transpose([[-1,1,1],[1,-1,1],[1,1,-1]])
                else:
                    lat_fam = 3
                    lat_type = 'rhombohedral'
                    niggli_num = 4
                    basis = np.transpose([[-1, 0, -1],[0, -1.32288, 0.5],
                                          [-0.548584, 0.774292, 1.04858]])

            elif np.allclose(2.0*abs(D+E+F),A+B,atol=eps):
                if np.allclose(D,E,atol=eps):
                    lat_fam = 4
                    lat_type = 'body centered tetragonal'
                    niggli_num = 6
                    basis = np.transpose([[-1, 1, 2], [1, 1.60788, -1.55394],
                                          [1.80278, -1.47253, 0.762655]])
                elif np.allclose(E,F,atol=eps):
                    lat_fam = 4
                    lat_type = 'body centered tetragonal'
                    niggli_num = 7
                    basis = np.transpose([[-1.95095, 1.41625, -0.433603], [1, -1, -2],
                                          [1.95095, 1.19163, 0.879663]])
                else:
                    lat_fam = 5
                    lat_type = 'body centered orthorhombic'
                    niggli_num = 8
                    basis = np.transpose([[1.41144, 0.0885622, -2],
                                          [-0.99868, 2.21232, 1.268178],
                                          [3.41012, -1.1237578, -1.268178]])
                    
    if np.allclose(A,B,atol=eps) and niggli_num is None:
        if positive:
            if np.allclose(D,E,atol=eps) and np.allclose(D,F,atol=eps) and np.allclose(A/2.,D,atol=eps):
                lat_fam = 3
                lat_type = 'rhombohedral'
                niggli_num = 9
                basis = np.transpose([[1,2,2],[2,1,2],[4,3,3]])
            elif np.allclose(D,E,atol=eps):
                lat_fam = 6
                lat_type = 'base centered monoclinic'
                niggli_num = 10
                basis = np.transpose([[-1.46391, 0, 1.96391], [1, 1, 1],
                                      [0, 2, 0]])
        else:
            if np.allclose(D,E,atol=eps) and np.allclose(D,F,atol=eps) and np.allclose(D,0,atol=eps):
                lat_fam = 4
                lat_type = 'simple tetragonal'
                niggli_num = 11
                basis = np.transpose([[1,0,0],[0,1,0],[0,0,2]])
            elif np.allclose(D,E,atol=eps):
                if np.allclose(D,0,atol=eps) and np.allclose(-A/2.,F,atol=eps):
                    lat_fam = 2
                    lat_type = 'hexagonal'
                    niggli_num = 12
                    basis = np.transpose([[1,0,0],[0.5,-0.8660254037844386,0],[0,0,2]])
                elif np.allclose(D,-A/2.,atol=eps) and np.allclose(F,0,atol=eps):
                    lat_fam = 4
                    lat_type = 'body centered tetragonal'
                    niggli_num = 15
                    basis = np.transpose([[-1,1,2],[1,-1,2],[1,1,-2]])
                elif np.allclose(D,0,atol=eps):
                    lat_fam = 5
                    lat_type = 'base centered orthorhombic'
                    niggli_num = 13
                    basis = np.transpose([[1, 1, 1], [1, -1, -1],
                                          [0, -1.73205, 1.73205]])
                elif np.allclose(2.0*abs(D+E+F),A+B,atol=eps):
                    lat_fam = 5
                    lat_type = 'face centered orthorhombic'
                    niggli_num = 16                    
                    basis = np.transpose([[ 1.04442 ,  1.43973 ,  1.68415 ],
                                          [ 0.779796, -1.1798  ,  1.      ],
                                          [ 1.779796, -0.1798  ,  0.      ]])
                else:
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 14
                    basis = np.transpose([[1,1,0],[0,2,0],[0.5,0,2]])
            elif np.allclose(2.0*abs(D+E+F),A+B,atol=eps):
                lat_fam = 6
                lat_type = 'base centered monoclinic'
                niggli_num = 17
                basis = np.transpose([[-0.05387, -0.61088, 2.51474], [1, 1, 1],
                                      [1.809568, -0.15957, 0]])

    if np.allclose(B,C,atol=eps) and niggli_num is None:
        if positive:
            if np.allclose(E,F,atol=eps):
                if np.allclose(A/4.,D,atol=eps) and np.allclose(A/2.,E,atol=eps):
                    lat_fam = 4
                    lat_type = 'body centered tetragonal'
                    niggli_num = 18
                    basis = np.transpose([[0, 0, 2], [1, -2, 1], [-2, -1, 1]])
                elif np.allclose(A/2.,E,atol=eps):
                    lat_fam = 5
                    lat_type = 'body centered orthorhombic'
                    niggli_num = 19
                    basis = np.transpose([[0.5, 1, 1.5], [0, 2, 0], [0, 0, 3]])
                else:
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 20
                    basis = np.transpose([[1, 1, 1], [1.70119, -1.45119, 1],
                                          [0.69779, -1.4322505, 3.23446]])
        else:
            if np.allclose(E,F,atol=eps):
                if np.allclose(D,0,atol=eps) and np.allclose(E,0,atol=eps):
                    lat_fam = 4
                    lat_type = 'simple tetragonal'
                    niggli_num = 21
                    basis = np.transpose([[0,0,0.5],[1,0,0],[0,1,0]])
                elif np.allclose(D,-B/2.,atol=eps) and np.allclose(E,0,atol=eps):
                    lat_fam = 2
                    lat_type = 'hexagonal'
                    niggli_num = 22
                    basis = np.transpose([[0, 0, -0.5], [1, 0, 0],
                                          [-0.5, 0.8660254037844386, 0]])
                elif np.allclose(E,0,atol=eps):
                    lat_fam = 5
                    lat_type = 'base centered orthorhombic'
                    niggli_num = 23
                    basis = np.transpose([[-0.3333333, -1.54116, 1.87449], [1, 1, 1],
                                          [2, -1, -1]])
                elif np.allclose(2.0*abs(D+E+F),A+B,atol=eps) and np.allclose(E,-A/3.,atol=eps):
                    lat_fam = 3
                    lat_type = 'rhombohedral'
                    niggli_num = 24
                    basis = np.transpose([[-1, 0, -1], [1.51184, 0, -0.845178],
                                          [-0.255922, -1.44338, 0.92259]])
                else:
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 25
                    basis = np.transpose([[1, 1, 1], [1.45119, -1.70119, -1],
                                          [0.28878, -3.26895, 0.48018]])

    if niggli_num is None:
        if positive:
            if np.allclose(E,F,atol=eps):
                if np.allclose(D,A/4.,atol=eps) and np.allclose(A/2.,E,atol=eps):
                    lat_fam = 5
                    lat_type = 'face centered orthorhombic'
                    niggli_num = 26
                    basis = np.transpose([[0, 1, 1.5], [0.5, 0, 1.5], [0, 0, 3]])
                elif np.allclose(A/2.,E,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 27
                    basis = np.transpose([[0.464824, -1.464824, -1.907413],
                                          [-1.618033, 0.618033, -1],
                                          [-1, -1, 0]])                    
            else:
                if np.allclose(E,A/2.,atol=eps) and np.allclose(F,2*D,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 28
                    basis = np.transpose([[-1.44896, 0.948958, -1], [-1, -1, 0],
                                          [0.342424, -1.342424, -2.02006]])                    
                elif np.allclose(F,A/2.,atol=eps) and np.allclose(E,2*D,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 29
                    basis = np.transpose([[-0.666125, 1.16613, 2.04852], [1, 1, 0],
                                          [1.61803, -0.618034, 1]])
                elif np.allclose(D,B/2.,atol=eps) and np.allclose(F,2*E,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 30
                    basis = np.transpose([[1, 1, 0], [1.61803, -0.618034, 1],
                                          [-0.0361373, 0.536137, 2.38982]])
                else:
                    lat_fam = 7
                    lat_type = 'triclinic'
                    niggli_num = 31
                    basis = None
        else:
            if np.allclose(E,F,atol=eps) and np.allclose(E,0,atol=eps):
                if np.allclose(D,0,atol=eps):
                    lat_fam = 5
                    lat_type = 'simple orthorhombic'
                    niggli_num = 32
                    basis = np.transpose([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
                elif np.allclose(D,-B/2.,atol=eps):
                    lat_fam = 5
                    lat_type = 'base centered orthorhombic'
                    niggli_num = 40
                    basis = np.transpose([[1, 1, 1], [1.61803, -0.618034, -1],
                                          [-1.05557, 1.99895, -0.943376]]) 
                else:
                    lat_fam = 6
                    lat_type = 'simple monoclinic'
                    niggli_num = 35
                    basis = np.transpose([[1,1,1],[1.61803,-0.618034,-1],
                                          [-0.668912,1.96676,-1.29785]])
                    
            elif np.allclose(D,F,atol=eps) and np.allclose(D,0,atol=eps):
                if np.allclose(E,-A/2.,atol=eps):
                    lat_fam = 5
                    lat_type = 'base centered orthorhombic'
                    niggli_num = 36
                    basis = np.transpose([[1, 1, 1], [1.41421, -1.41421, 0],
                                          [-1.43541, -1.43541, 1.37083]]) 
                else:
                    lat_fam = 6
                    lat_type = 'simple monoclinic'
                    niggli_num = 33
                    basis = np.transpose([[2,0,0],[0,2,0],[0.5,0,2]])
            elif np.allclose(D,E,atol=eps) and np.allclose(D,0,atol=eps):
                if np.allclose(F,-A/2.,atol=eps):
                    lat_fam = 5
                    lat_type = 'base centered orthorhombic'
                    niggli_num = 38
                    basis = np.transpose([[0.5, 1, 0], [0.5, -1, 0], [0, 0, 3]]) 
                else:
                    lat_fam = 6
                    lat_type = 'simple monoclinic'
                    niggli_num = 34
                    basis = np.transpose([[1,1,1],[1.22474487,-1.22474487,-2],
                                          [-0.16598509,-1.64308297,1.80906806]])
            else:
                if np.allclose(-B/2.,D,atol=eps) and np.allclose(-A/2.,E,atol=eps) and np.allclose(F,0,atol=eps):
                    lat_fam = 5
                    lat_type = 'body centered orthorhombic'
                    niggli_num = 42
                    basis = np.transpose([[-1.53633, 1.36706, -1.33073], [1, 1, 1],
                                          [1.61803, -0.61803, -1]]) 
                elif np.allclose(-B/2.,D,atol=eps) and np.allclose(F,0,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 41
                    basis = np.transpose([[-1.85397, -0.854143, 1.35397],[1, 0, 1],
                                          [1, -1.41421, -1]]) 
                elif np.allclose(E,-A/2.,atol=eps) and np.allclose(F,0,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 37
                    basis = np.transpose([[-1.79092,-1.47209,0.790922],[1.0,-1.41421,-1.0],
                                          [1.0,0.0,1.0]]) 
                elif np.allclose(E,0,atol=eps) and np.allclose(F,-A/2.,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 39
                    basis = np.transpose([[0, -1.73205,-1],[-1.66542, -0.672857, 1.66542],
                                           [1,0,1]]) 
                elif np.allclose(2.0*abs(D+E+F),A+B,atol=eps) and np.allclose(abs(2.0*D+F),B,atol=eps):
                    lat_fam = 6
                    lat_type = 'body centered monoclinc'
                    niggli_num = 43
                    basis = np.transpose([[-0.39716, -0.34718, 2.49434],
                                          [2.64194, -0.14194, 0],
                                          [-1.39716, -1.34718, 1.49434]]) 
                else:
                    lat_fam = 7
                    lat_type = 'triclinic'
                    niggli_num = 44
                    basis = None

    return lat_type, niggli_num, lat_fam, basis 
