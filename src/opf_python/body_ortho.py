def body_ortho_srHNFs(n):
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
    n. Assuming A = [[ 1.  , 1.  , 2.  ], [ 1.41144 , 0.0885622, -2.
    ], [-1.99868 , 1.21232 , -0.731822 ]].

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
        
        if f%a==0 and f%c==0:
            if a==c and c==f:
                es = [0]
                ds = [0]
                bs = [0]
            elif a==f and a!=c:
                bs = [0]
                ds = [0]
                es = range(0,f,c)
            elif c==f and a!=c:
                es = [0]
                bs = range(c)
                ds = range(f)
            elif a*c*f == f:
                bs = [0]
                ds = find_es(f)
                es = ds
            else:
                es = range(0,f,c)
                bs = range(c)
                ds = range(f)
                
            for e in es:
                a23 = -c+e
                g21 = c-e*e/float(c)
                if a23%a==0 and g21%f==0:
                    for b in bs:
                        b23 = e-b*a23/float(a)
                        b33 = f-b*f/float(a)
                        if b33%c==0 and b23%c==0:
                            for d in ds:
                                a13 = -b+d
                                b11 = -a+b+d
                                b13= -a+d-b*a13/float(a)
                                g11 = b11-e*b11/float(c)
                                g13 = d-d*a13/float(a)-e*b13/float(c)
                                g23 = e-d*a23/float(a)-e*b23/float(c)
                                if a13%a==0 and b11%c==0 and b13%c==0 and g11%f==0 and g13%f==0 and g23%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    spHNFs.append(HNF)
                                
    return spHNFs
                                    
def find_es(f):
    """Finds the allowed values of e (and d) for the niggli 8 cells.

    Args:
        f (int): The value of the lower right corner of the HNF.
    
    Returns:
        es (list): The allowed values of e and d.
    """
    import numpy as np
    smallest = (1-(f-1)**2)/f
    ops = range(smallest,1)
    es = [np.sqrt(-(o*f)+1) for o in ops if np.sqrt(-(o*f)+1)%1==0]
    if len(es) > 0:
        return sorted(es)
    else: 
        return [0]

def body_ortho_19(n):
    """Finds the symmetry preserving HNFs for the body centered
    orthorhombic lattices with niggli setting 19 a determinant of
    n. Assuming A = [[-1. , 0. , 0. ], [-0.5, -1. , 1.5], [-0.5, 1. ,
    1.5]].

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
        
        if f%a==0:# and f%c==0:
            es = range(f)
            ds = range(f)
            bs = range(c)
            for b in bs:
                b32 = b*f/float(a)
                b31 = f+b32
                if b32%c==0 and b31%c==0:
                    for e in es:
                        a21 = -c-e
                        b21 = -b*a21/float(a)+e
                        b22 = -b*a21/float(a)
                        if a21%a==0 and b21%c==0 and b22%c==0:
                            for d in ds:
                                a11 = -b-d
                                b11 = b-b*a11/float(a)+d
                                b12 = 2*b-b*a11/float(a)
                                g11 = b+d-a11*d/float(a)-b11*e/float(c)
                                g12 = 2*d -a11*d/float(a) -b12*e/float(c)
                                g21 = c-d*a21/float(a)-e*b21/float(c)
                                g22 = -d*a21/float(a)+e*b22/float(c)
                                if a11%a==0 and b11%c==0 and b12%c==0 and g11%f==0 and g12%f==0 and g21%f==0 and g22%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    spHNFs.append(HNF)
                                
    return spHNFs                
