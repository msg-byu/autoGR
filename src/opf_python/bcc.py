def bcc_srHNFs(n):
    """Finds the symmetry preserving HNFs for the base centered cubic lattices
    with a determinant of n. Assuming A = [[-1,1,1],[1,-1,1],[1,1,-1]]


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
                for b in range(c):
                    #find e values alpha2 condition
                    if (b*f)%(a*c)==0:
                        e = 0
                        while (e<f):
                            #beta2 condition
                            g22 = -c +e*e/float(c)
                            beta21 = b*e/float(a)
                            if e%c==0 and g22%f ==0 and beta21%c==0:
                                d = 0
                                #alpha1 condition
                                while (d<f):
                                    beta11 = -a+b-b*d/float(a)
                                    beta12 = b-d
                                    g11 = -a+b-d*(d/float(a)-1)-e*beta11/float(c)
                                    g12 = -b+d-beta12*e/float(c)
                                    g21 = c-e*e*b/float(a*c) -d*e/float(a)
                                    if beta11%c==0 and beta12%c==0 and g11%f==0 and g12%f==0 and g21%f==0:
                                        HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                        srHNFs.append(HNF)
                                
                                
                                    d += a
                            e += a
    return srHNFs
