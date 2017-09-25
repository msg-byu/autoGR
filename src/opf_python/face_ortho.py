def face_ortho_srHNFs(n):
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

def face_ortho_26(n):
    """Finds the symmetry preserving HNFs for the face centered
    orthorhombic lattices with a niggli cell of 26 and a determinant
    of n. Assuming A = [[ 1. , 0. , 0. ], [ 0.5, 1. , 0. ], [ 0.5,
    0. , 1.5]].

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

        if f%a==0 and c%a==0:
            bs = range(0,c,a)#[0]
            ds = range(0,f,a)
            es = range(0,f,a)

            for b in bs:
                b12 = 2*b+b*b/float(a)
                if b*f%(c*a)==0 and b12%c==0:
                    for e in es:
                        a21 = -c-e
                        b21 = b*a21/float(a)
                        b23 = b*e/float(a)
                        if a21%a==0 and b21%c==0 and b23%c==0:
                            for d in ds:
                                a11 = -b-d
                                b11 = 2*b-b*a11/float(a)
                                b13 = b*d/float(a)
                                g11 = 2*d-d*a11/float(a)-b11*e/float(c)
                                g12 = b13-b12*e/float(c)
                                g13 = 2*d+d*d/float(a)-b13*e/float(c)
                                g21 = -d*a21/float(a)+e*b21/float(c)
                                g22 = c*d/float(a)-2*e-b*e/float(a)
                                g23 = 2*e+d*e/float(a)-e*b23/float(c)
                                if a11%a==0 and b11%c==0 and b13%c==0 and g11%f==0 and g12%f==0 and g13%f==0 and g21%f==0 and g22%f==0 and g23%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    spHNFs.append(HNF)
                                
    return spHNFs
                                    
                                
            
