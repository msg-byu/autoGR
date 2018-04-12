def base_mono_20_25(n):
    """Finds the symmetry preserving HNFs for the base centered monoclinic
    lattices with a determinant of n. Assuming for basis 20 A = [[ 1.
    , 1.  , 1.  ], [ 1.70119 , -1.45119 , 1.  ], [ 0.69779 ,
    -1.4322505, 3.23446 ]], for basis 25 A = [[ 1.  , 1.  , 1.  ], [
    1.45119, -1.70119, -1.  ], [ 0.28878, -3.26895, 0.48018]].

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        spHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    spHNFs = []
    
    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]

        if c%2==0:
            bs = [0,c//2.]
        else:
            bs = [0]

        for e in range(f):
            g22 = -c-2*e
            if g22%f==0:
                for b in bs:
                    g12 = -b-2*b*e/float(c)
                    if g12%f==0:
                        for d in range(f):
                            HNF = [[a,0,0],[b,c,0],[d,e,f]]
                            spHNFs.append(HNF)
                                
    return spHNFs

def base_mono_28(n):
    """Finds the symmetry preserving HNFs for the base centered monoclinic
    lattices with a determinant of n. Assuming A = [[-1.44896 ,
    0.948958, -1.  ], [-1.  , -1.  , 0.  ], [ 0.342424, -1.342424,
    -2.02006 ]].

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        spHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    spHNFs = []
    
    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]

        if f%c==0:
            for e in range(0,f,c):
                g21 = 2*e+e*e/float(c)
                if g21%f==0:
                    for d in range(0,f,c):
                        g11 = 2*d+d*e/float(c)
                        if g11%f==0:
                            for b in range(c):
                                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                spHNFs.append(HNF)
                                
    return spHNFs

def base_mono_29_30(n):
    """Finds the symmetry preserving HNFs for the base centered monoclinic
    lattices with a determinant of n. Assuming for basis 29 A =
    [[-0.666125, 1.16613 , 2.04852 ], [ 1.  , 1.  , 0.  ], [ 1.61803 ,
    -0.618034, 1.  ]], for basis 30 A = [[ 1.  , 1.  , 0.  ], [
    1.61803 , -0.618034 , 1.  ], [-0.0361373, 0.536137 , 2.38982 ]].

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        spHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    spHNFs = []
    
    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]

        if f%c==0:
            for e in range(0,f,c):
                g21 = 2*e+e*e/float(c)
                if g21%f==0:
                    for d in range(0,f,c):
                        g11 = 2*d+d*e/float(c)
                        if g11%f==0:
                            for b in range(c):
                                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                spHNFs.append(HNF)
                                
    return spHNFs
                                                          
def base_mono_10_14_17_27_37_39_41(n):
    """Finds the symmetry preserving HNFs for the base centered monoclinic
    lattices with a determinant of n. Assuming:
	for basis 10 A = [[1, -1, 1],[-1.46391, 0, 1.96391],[0, 2, 0]],
	for basis 14 A = [[-1,1,0],[0.5,0,2],[0,-2,0]],
	for basis 17 A = [[-1.05387,-1.61088,1.51474],[-0.244302,-2.77045,0.51474],[1.809568,-0.15957,0.0]]
	for basis 27 A = [[0.0,-1.73205,-1.0],[-1.66542,-0.672857,1.66542],[1.0,0.0,1.0]],
	for basis 37 A = [[-1.79092,-1.47209,0.790922],[1.0,-1.41421,-1.0],[1.0,0.0,1.0]],
	for basis 39 A = [[0, -1.73205,-1],[-1.66542, -0.672857, 1.66542], [1,0,1]],
	for basis 41 A = [[-1.85397, -0.854143, 1.35397],[1, 0, 1],[1, -1.41421, -1]].

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        spHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals
    
    diags = get_HNF_diagonals(n)

    spHNFs = []

    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]
        #gamma 2
        if f%2==0:
            es = [0,f//2.]
        else:
            es = [0]

        for e in es:
            for d in range(f):
                g11 = -a+2*d
                if g11%f==0:
                    for b in range(c):
                        HNF = [[a,0,0],[b,c,0],[d,e,f]]
                        spHNFs.append(HNF)

    return spHNFs 
                      
def base_mono_43(n):
    """Finds the symmetry preserving HNFs for the base centered monoclinic
    lattices with a determinant of n. Assuming A = [[-0.39716,
    -0.34718, 2.49434], [ 2.64194, -0.14194, 0.  ], [-1.39716,
    -1.34718, 1.49434]].

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        spHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    spHNFs = []
    
    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]

        if f%c==0:
            for e in range(0,f,c):
                g22 = 2*e-e*e/float(c)
                if g22%f==0:
                    for d in range(f):
                        b12 = a+d
                        g12 = 2*a+2*d-b12*e/float(c)
                        if b12%c==0 and g12%f==0:
                            for b in range(c):
                                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                spHNFs.append(HNF)
                                
    return spHNFs
