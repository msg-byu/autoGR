def body_ortho_19(n):
    """Finds the symmetry preserving HNFs for the body centered
    orthorhombic lattices with a determinant of n. Assuming A =
    [[0.5,1,1.5],[0,2,0],[0,0,3]] (3rd basis choince in
    notes/body_centered_ortho).

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
        
        #beta1 condition
        if c%2==0:
            bs = [0,c/2]
        else:
            bs = [0]
        #gamma1 and gamma2 conditions
        if f%2==0:
            es = [0,f/2]
        else:
            es = [0]

        for b in range(c):
            #beta1 condition
            beta13 = a+2*b
            if beta13%c==0:
                for e in es:
                    #gamma1 condition
                    gamma13 = e*beta13/float(c)
                    if gamma13%f==0:
                        for d in range(f):
                            #gamma1 condition
                            gamma12 = a+2*d
                            if gamma12%f==0:
                                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                spHNFs.append(HNF)
                                
    return spHNFs

def body_ortho_8(n):
    """Finds the symmetry preserving HNFs for the body centered
    orthorhombic lattices with niggli setting 8 a determinant of
    n. Assuming A = [[ 1.41144 , 0.0885622, -2.  ], [-0.99868 ,
    2.21232 , 1.268178 ], [ 3.41012 , -1.1237578, -1.268178 ]].

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

        if 2*f%c==0:
            if c%2==0:
                bs = [0,c//2.]
            else:
                bs = [0]

            for b in bs:
                for e in range(f):
                    b21 = 2*e
                    g21 = -b21+b21*e/float(c)
                    if b21%c==0 and g21%f==0:
                        for d in range(f):
                            b11 = -a+2*b-2*d
                            g11 = -b11*e/float(c)
                            g13 = a+2*d-b21*b/float(c)
                            if b11%c==0 and g11%f==0 and g13%f==0:
                                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                spHNFs.append(HNF)
                                
    return spHNFs

def body_ortho_42(n):
    """Finds the symmetry preserving HNFs for the body centered
    orthorhombic lattices with niggli setting 42 a determinant of
    n. Assuming A = [[-1.53633, 1.36706, -1.33073], [ 1.  , 1.  , 1.
    ], [ 1.61803, -0.61803, -1.  ]].

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

        if f%2==0:
            es = [0,f//2]
        else:
            es = [0]

        for e in es:
            for b in range(c):
                b11 = -a+2*b
                g12 = -b11*e/float(c)
                if b11%c==0 and g12%f==0:
                    for d in range(f):
                        g11 = -a+2*d-b11*e/float(c)
                        g13 = -a+2*d
                        if g11%f==0 and g13%f==0:
                            HNF = [[a,0,0],[b,c,0],[d,e,f]]
                            spHNFs.append(HNF)
                                
    return spHNFs
