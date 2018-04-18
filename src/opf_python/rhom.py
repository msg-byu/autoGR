def rhom_9_24(n):
    """Finds the symmetry preserving HNFs for the trigonal lattices
    with a determinant of n. Assuming A = [[1,2,2],[2,1,2],[4,3,3]] for basis 9
	A = [[-0.255922,-1.44338,0.92259],[1.51184,0,-0.845178],[1.255922,1.44338,0.07741]] 
	for basis 24.

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

        #alpha3 condition
        if f%a==0:
            #beta1 condition
            if c%2==0:
                bs = [0,c/2]
            else:
                bs = [0]
                    
            for b in bs:
                #beta3 condition
                beta13 = f+b*f/float(a)
                if beta13%c==0:                    
                    #find e values alpha2 condition
                    e = 0
                    while (e<f):
                        beta22 = e+b*e/float(a)
                        g21= c+2*e
                        g11 = b+2*b*e/float(c)
                        if beta22%c==0 and g21%f==0 and g11%f==0:
                            #find d from alpha1 condition
                            d = 0
                            while (d<f):
                                #beta2 condition
                                beta12 = -a+b+d+d*b/float(a)
                                g12 = -b-d+d*d/float(a)-e*beta12/float(c)
                                g22 = -c-2*e+d*e/float(a)-e*beta22/float(c)
                                if beta12%c==0 and g12%f==0 and g22%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    spHNFs.append(HNF)
                                
                                
                                d += a
                        e += a
    return spHNFs

def rhom_4_2(n):
    """Finds the symmetry preserving HNFs for the rhombohedral lattice
    that matches niggli conditions for niggli cells number 2 and 4
    with a determinant of n. Assuming A =
    [[-1.11652,-0.610985,0.616515],[0.0,-1.32288,-0.5],[1.0,1.32288,1.5]]
    for basis 2 and A =
    [[-0.548584,0.774292,1.04858],[0.0,-1.32288,0.5],[1.0,1.32288,0.5]]
    for basis 4.

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

        # conditions alpha 3 and betta3
        if f%a == 0:
            if c%2==0:
                bs = [0,c/2]
            else:
                bs = [0]

            for b in bs:
                beta32 = -f+b*f/float(a)
                if beta32%c==0:
                    for e in range(f):
                        beta22 = -e+b*e/float(a)
                        gamma11 = b-2*b*e/float(c)
                        gamma21 = c-2*e
                        if e%a==0 and beta22%c==0 and gamma11%f==0 and gamma21%f==0:
                            for d in range(f):
                                beta12 = -a+b-d+b*d/float(a)
                                gamma12 = -a+d*d/float(a)-e*beta12/float(c)
                                gamma22 = -e+d*e/float(a)-e*beta22/float(c)
                                if d%a==0 and beta12%c==0 and gamma12%f==0 and gamma22%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    spHNFs.append(HNF)
                                
    return spHNFs
                  
