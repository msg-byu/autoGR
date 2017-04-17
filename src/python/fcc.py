def fcc_srHNFs(n):
    """Finds the symmetry preserving HNFs for the face centered cubic lattices
    with a determinant of n. Assuming A = [[0,1,1],[1,0,1],[1,1,0]].

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
                #find b values beta3 condition
                b = 0
                while (b<c):
                    #find e values beta2 condition
                    e = 0
                    while (e<f):
                        #alpha2 condition
                        alpha21 = c+e
                        beta21 = b*(e+c)/float(a)
                        g22 = c+e+e*e/float(c)
                        if alpha21%a==0 and beta21%c==0 and g22%f==0:
                            d = 0
                            #beta1 condition
                            while (d<f):
                                alpha11 = b+d
                                beta12 = b-d
                                beta11 = 2*b+b*(b+d)/float(a)
                                g11 = d*(2-alpha11/float(a))-e*beta11/float(c)
                                g12 = a+b+2*d-e*beta12/float(c)
                                g21 = -d*alpha21/float(a)-e*beta21/float(c)
                                if beta11%c==0 and beta12%c == 0 and alpha11%a==0 and g11%f==0 and g12%f==0 and g21%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    srHNFs.append(HNF)
                                
                                
                                d += a
                        e += c
                    b += a
    return srHNFs
