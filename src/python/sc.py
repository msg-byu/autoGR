def sc_srHNFs(n):
    """Finds the symmetry preserving HNFs for the simple cubic lattices
    with a determinant of n. Assuming A = [[1,0,0],[0,1,0],[0,0,1]].

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

        #alpha3 condition
        if f%a==0:
            #beta3 condition
            if f%c==0:
                #find b values beta2 condition
                b = 0
                while (b<c):
                    #find e values alpha2 condition
                    e = 0
                    while (e<f):
                        #beta2 condition
                        g21 = -c+e*e/float(c)
                        if e%c==0 and g21%f ==0:
                            d = 0
                            #alpha1 condition
                            while (d<f):
                                beta11 = b-d
                                beta12 = -a-b*d/float(a)
                                g11 = -b-beta11*e/float(c)+d
                                g12 = b-beta12*e/float(c)-d*d/float(a)
                                g22 = c-d*e/float(a) +b*e*e/float(a*c)
                                if beta11%c==0 and beta12%c==0 and g11%f==0 and g12%f==0 and g22%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    srHNFs.append(HNF)
                                
                                
                                d += a
                        e += a
                    b += a
    return srHNFs
