def base_ortho_38_13(n):
    """Finds the symmetry preserving HNFs for the base centered
    orthorhombic lattices with a determinant of n. Assuming A =
    [[0.5,1,0],[0.5,-1,0],[0,0,3]] (1st basis choince in
    notes/base_ortho) for case 38 and A = [[ 1.  , 1.  , 1.  ], [ 1.
    , -1.  , -1.  ], [ 0.  , -1.73205, 1.73205]] for case 13.

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
        #alpha2 condition
        if c%a==0:
            #gamma1 and gamma2 conditions
            if f%2==0:
                ed_vals = [0,f/2]
            else:
                ed_vals = [0]
                
            b=0
            #alhpa1 condition
            while (b<c):
                #beta1 condition
                beta13=-a+b*b/float(a)
                if beta13%c==0:
                    for e in ed_vals:
                        for d in ed_vals:
                            #gamma1 and gamma2 conditions
                            g13 = -d+b*d/float(a) -e*beta13/float(c)
                            g23 = c*d/float(a)-e-b*e/float(a)
                            if g13%f==0 and g23%f==0:
                                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                spHNFs.append(HNF)
                b += a
                                
    return spHNFs

def base_ortho_23(n):
    """Finds the symmetry preserving HNFs for the base centered
    orthorhombic lattices with a determinant of n. Assuming A =
    [[-0.3333333, -1.54116 , 1.87449 ], [ 1.  , 1.  , 1.  ], [ 2.  ,
    -1.  , -1.  ]].

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

        if f%a==0:
            if c%2==0:
                bs = [0,c//2.]
            else:
                bs = [0]
            if f%2==0:
                es = [0,f//2.]
            else:
                es = [0]

            for b in bs:
                if b*f%(a*c)==0:
                    for e in es:
                        if b*e%(a*c)==0 and 2*b*e%(c*f)==0 and e%a==0:
                            for d in range(0,f,a):
                                b13 = -b+b*d/float(a)
                                g13 = -a+d*d/float(a)-b13*e/float(c)
                                g23 = e+d*e/float(a)-b*e*e/float(a*c)
                                if b13%c==0 and g13%f==0 and g23%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    spHNFs.append(HNF)
                                
    return spHNFs
                                                                
def base_ortho_40(n):
    """Finds the symmetry preserving HNFs for the base centered
    orthorhombic lattices with a determinant of n. Assuming A = [[ 1.
    , 1.  , 1.  ], [ 1.61803 , -0.618034, -1.  ], [-1.05557 , 1.99895
    , -0.943376]].

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
                    for d in range(0,f,c):
                        g12 = 2*d-d*e/float(c)
                        if g12%f==0:
                            for b in range(c):
                                b13 = 2*b-d
                                g13 = b13*e/float(c)
                                if b13%c==0 and g13%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    spHNFs.append(HNF)
                                
    return spHNFs
                                
def base_ortho_36(n):
    """Finds the symmetry preserving HNFs for the base centered
    orthorhombic lattices with a determinant of n. Assuming A = [[1,
    1, 1], [1.41421, -1.41421, 0], [-1.43541, -1.43541, 1.37083]].

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

        if f%a==0:
            if f%2==0:
                es = [0,f//2.]
            else:
                es = [0]
            if c%2==0:
                bs = [0,c//2.]
            else:
                bs = [0]

            for b in bs:
                b32 = -b*f/float(a)
                if b32%c==0:
                    for e in es:
                        b22 = -b*e/float(a)
                        g13 = -2*b*e/float(c)
                        if b22%c==0 and g13%f==0 and e%a==0:
                            for d in range(0,f,a):
                                b12 = -b*d/float(a)
                                g12 = 2*d-d*d/float(a)-b12*e/float(c)
                                g22 = 2*e-d*e/float(a)-b22*e/float(c)
                                if b12%c==0 and g12%f==0 and g22%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    spHNFs.append(HNF)
                                
    return spHNFs                                    
