def hex_srHNFs(n):
    """Finds the symmetry preserving HNFs for the hexagonal lattices
    with a determinant of n. Assuming A = [[1,0,0],[0.5,-0.8660254037844386,0],[0,0,2]].

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        srHNFs (list of lists): The symmetry preserving HNFs.
    """

    from universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    srHNFs = []
    
    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]

        #alpha2 condition
        if c%a==0:
            b=0
            #alpha1 condition
            while (b<c):
                beta13 = a+2*b
                beta11 = 2*b+b*b/float(a)
                if beta13%c==0 and beta11%c==0:
                    #gamma2 condition
                    if f%2==0:
                        es = [0,f/2]
                    else:
                        es = [0]

                    for e in es:
                        gamma13 = (a+2*b)*e/float(c)
                        if gamma13%f==0:
                            for d in range(f):
                                gamma11 = b*d/float(a)-e*beta11/float(c)
                                gamma12 = 2*d+b*d/float(a) - e*beta11/float(c)
                                gamma21 = c*d/float(a) -2*e -b*e/float(a)
                                gamma22 = (c*d-b*e)/float(a)
                                if gamma11%f==0 and gamma12%f==0 and gamma21%f==0 and gamma22%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    srHNFs.append(HNF)
                b += a
                                
    return srHNFs
