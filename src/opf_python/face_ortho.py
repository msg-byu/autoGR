def face_ortho_26(n):
    """Finds the symmetry preserving HNFs for the face centered
    orthorhombic lattices with a determinant of n. Assuming A =
    [[0,1,1.5],[0.5,0,1.5],[0,0,3]] (10th basis choince in
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

        #beta1 candition
        if c%2==0:
            bs = [0,c/2]
        else:
            bs = [0]

        for b in bs:
            for e in range(f):
                #gamma2 condition and gamma1 condition
                gamma23 = c+2*e
                gamma12 = b+2*b*e/float(c)
                if gamma23%f==0 and gamma12%f==0:
                    for d in range(f):
                        #gamma1 condition
                        gamma13 = a+b+2*d
                        if gamma13%f==0:
                            HNF = [[a,0,0],[b,c,0],[d,e,f]]
                            spHNFs.append(HNF)
                                
    return spHNFs

def face_ortho_16(n):
    """Finds the symmetry preserving HNFs for the face centered
    orthorhombic lattices with a niggli cell of 26 and a determinant
    of n. Assuming A = [[1.04442, 1.43973, 1.68415], [0.779796, -1.1789, 1.0
    ], [1.779796, 0.1789, 0]].

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        spHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    spHNFs=[]

    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]
        
        #beta 1 condition
        if c%2==0:
            bs = [c/2,0]
        else:
            bs = [0]
            
        #gamma 1 condition
        for b in bs:
            for e in range(f):
                g11=-b-2*b*e/float(c)
                g21=c+2*e
                if g11%f==0 and g21%f==0:
                    for d in range(f):
                        g12=a+2*d-(2*b*e)/float(c)
                        if g12%f==0:
                            HNF = [[a,0,0],[b,c,0],[d,e,f]]
                            spHNFs.append(HNF)

    return spHNFs
