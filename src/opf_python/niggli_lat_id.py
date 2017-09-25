"""The methods contained in this module are used to identify a lattice
using the niggli reduced cell."""

import numpy as np

def niggli_id(lattice,eps_=None):
    """Returns the lattice type, niggli cell number, and crystal family
    number, .i.e, cubic -> 1, hexagonal -> 2, rhombehedral -> 3,
    tetragonal -> 4, orthorhombic -> 5, monoclinic -> 6, and triclinic
    -> 7. 

    Args:
        lattice (numpy array): The lattice vectors as columns of a matrix.
        eps_ (float, optional): The floating point tollerance for comparisons.

    Returns:
        lat_type, niggli_num, lat_fam (str, int, int): The lattice type, the 
            niggli case number from the 'International Tables of Crystalography',
            the lattice family type. 
    """

    from opf_python.pyniggli import reduced_cell

    if not isinstance(lattice,np.ndarray):
        lattice = np.array(lattice)

    if eps_ is None:
        eps = 1E-5
    else:
        eps = eps_
        
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

    positive = False
    if D > 0 and E > 0 and F > 0:
        positive = True

    niggli_num = None
    if np.allclose(A,B,atol=eps) and np.allclose(B,C,atol=eps):
        if positive:
            if np.allclose(D,E,atol=eps) and np.allclose(D,F,atol=eps):
                if np.allclose(A/2., D,atol=eps) and np.allclose(A/2., E,atol=eps) and np.allclose(A/2., F,atol=eps):
                    lat_fam = 1
                    lat_type = 'face centered cubic'
                    niggli_num = 1
                else:
                    lat_fam = 3
                    lat_type = 'rhombohedral'
                    niggli_num = 2
        else:
            if np.allclose(D,E,atol=eps) and np.allclose(D,F,atol=eps):
                if np.allclose(0,D,atol=eps):
                    lat_fam = 1
                    lat_type = 'simple cubic'
                    niggli_num = 3
                elif np.allclose(-A/3., D,atol=eps):
                    lat_fam = 1
                    lat_type = 'body centered cubic'
                    niggli_num = 5
                else:
                    lat_fam = 3
                    lat_type = 'rhombohedral'
                    niggli_num = 4

            elif np.allclose(2.0*abs(D+E+F),A+B,atol=eps):
                if np.allclose(D,E,atol=eps):
                    lat_fam = 4
                    lat_type = 'body centered tetragonal'
                    niggli_num = 6
                elif np.allclose(E,F,atol=eps):
                    lat_fam = 4
                    lat_type = 'body centered tetragonal'
                    niggli_num = 7
                else:
                    lat_fam = 5
                    lat_type = 'body centered orthorhombic'
                    niggli_num = 8
                    
    if np.allclose(A,B,atol=eps) and niggli_num is None:
        if positive:
            if np.allclose(D,E,atol=eps) and np.allclose(D,F,atol=eps) and np.allclose(A/2.,D,atol=eps):
                lat_fam = 3
                lat_type = 'rhombohedral'
                niggli_num = 9
            elif np.allclose(D,E,atol=eps):
                lat_fam = 6
                lat_type = 'bace centered monoclinic'
                niggli_num = 10
        else:
            if np.allclose(D,E,atol=eps) and np.allclose(D,F,atol=eps) and np.allclose(D,0,atol=eps):
                lat_fam = 4
                lat_type = 'simple tetragonal'
                niggli_num = 11
            elif np.allclose(D,E,atol=eps):
                if np.allclose(D,0,atol=eps) and np.allclose(-A/2.,F,atol=eps):
                    lat_fam = 2
                    lat_type = 'hexagonal'
                    niggli_num = 12
                elif np.allclose(D,-A/2.,atol=eps) and np.allclose(F,0,atol=eps):
                    lat_fam = 4
                    lat_type = 'body centered tetragonal'
                    niggli_num = 15
                elif np.allclose(D,0,atol=eps):
                    lat_fam = 5
                    lat_type = 'base centered orthorhombic'
                    niggli_num = 13
                elif np.allclose(2.0*abs(D+E+F),A+B,atol=eps):
                    lat_fam = 5
                    lat_type = 'face centered orthorhombic'
                    niggli_num = 16                    
                else:
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 14
            elif np.allclose(2.0*abs(D+E+F),A+B,atol=eps):
                lat_fam = 6
                lat_type = 'bace centered monoclinic'
                niggli_num = 17

    if np.allclose(B,C,atol=eps) and niggli_num is None:
        if positive:
            if np.allclose(E,F,atol=eps):
                if np.allclose(A/4.,D,atol=eps) and np.allclose(A/2.,E,atol=eps):
                    lat_fam = 4
                    lat_type = 'body centered tetragonal'
                    niggli_num = 18
                elif np.allclose(A/2.,E,atol=eps):
                    lat_fam = 5
                    lat_type = 'body centered orthorhombic'
                    niggli_num = 19
                else:
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 20
        else:
            if np.allclose(E,F,atol=eps):
                if np.allclose(D,0,atol=eps) and np.allclose(E,0,atol=eps):
                    lat_fam = 4
                    lat_type = 'simple tetragonal'
                    niggli_num = 21
                elif np.allclose(D,-B/2.,atol=eps) and np.allclose(E,0,atol=eps):
                    lat_fam = 4
                    lat_type = 'hexagonal'
                    niggli_num = 22
                elif np.allclose(E,0,atol=eps):
                    lat_fam = 5
                    lat_type = 'base centered orthorhombic'
                    niggli_num = 23
                elif np.allclose(2.0*abs(D+E+F),A+B,atol=eps) and np.allclose(E,-A/3.,atol=eps):
                    lat_fam = 3
                    lat_type = 'rhombohedral'
                    niggli_num = 24
                else:
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 25

    if niggli_num is None:
        if positive:
            if np.allclose(E,F,atol=eps):
                if np.allclose(D,A/4.,atol=eps) and np.allclose(A/2.,E,atol=eps):
                    lat_fam = 5
                    lat_type = 'face centered orthorhombic'
                    niggli_num = 26
                elif np.allclose(A/2.,E,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 27
            else:
                if np.allclose(E,A/2.,atol=eps) and np.allclose(F,2*D,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 28
                elif np.allclose(F,A/2.,atol=eps) and np.allclose(E,2*D,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 29
                elif np.allclose(D,B/2.,atol=eps) and np.allclose(F,2*E,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 30
                else:
                    lat_fam = 7
                    lat_type = 'triclinic'
                    niggli_num = 31
        else:
            if np.allclose(E,F,atol=eps) and np.allclose(E,0,atol=eps):
                if np.allclose(D,0,atol=eps):
                    lat_fam = 5
                    lat_type = 'simple orthorhombic'
                    niggli_num = 32
                elif np.allclose(D,-B/2.,atol=eps):
                    lat_fam = 5
                    lat_type = 'base centered orthorhombic'
                    niggli_num = 40
                else:
                    lat_fam = 6
                    lat_type = 'simple monoclinic'
                    niggli_num = 35

            elif np.allclose(D,F,atol=eps) and np.allclose(D,0,atol=eps):
                if np.allclose(E,-A/2.,atol=eps):
                    lat_fam = 5
                    lat_type = 'base centered orthorhombic'
                    niggli_num = 36
                else:
                    lat_fam = 6
                    lat_type = 'simple monoclinic'
                    niggli_num = 33
            elif np.allclose(D,E,atol=eps) and np.allclose(D,0,atol=eps):
                if np.allclose(F,-A/2.,atol=eps):
                    lat_fam = 5
                    lat_type = 'base centered orthorhombic'
                    niggli_num = 38
                else:
                    lat_fam = 6
                    lat_type = 'simple monoclinic'
                    niggli_num = 34
            else:
                if np.allclose(-B/2.,D,atol=eps) and np.allclose(-A/2.,E,atol=eps) and np.allclose(F,0,atol=eps):
                    lat_fam = 5
                    lat_type = 'body centered orthorhombic'
                    niggli_num = 42
                elif np.allclose(-B/2.,D,atol=eps) and np.allclose(F,0,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 41
                elif np.allclose(E,-A/2.,atol=eps) and np.allclose(F,0,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 37
                elif np.allclose(E,0,atol=eps) and np.allclose(F,-A/2.,atol=eps):
                    lat_fam = 6
                    lat_type = 'base centered monoclinic'
                    niggli_num = 39
                elif np.allclose(2.0*abs(D+E+F),A+B,atol=eps) and np.allclose(abs(2.0*D+F),B,atol=eps):
                    lat_fam = 6
                    lat_type = 'body centered monoclinc'
                    niggli_num = 43
                else:
                    lat_fam = 7
                    lat_type = 'triclinic'
                    niggli_num = 44

    return lat_type, niggli_num, lat_fam
