def body_tet_srHNFs(n):
    """Finds the symmetry preserving HNFs for the body centered tetragonal
    lattices with a determinant of n. Assuming A =
    [[-1,1,2],[1,-1,2],[1,1,-2]] (first basis choince in notes/body_centered_tet.nb).

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
            if f%c ==0:
                for b in range(c):
                    #beta3 condition
                    if (b*f)%(a*c) == 0:
                        e = 0
                        #alpha2 condition
                        while(e<f):
                            #beta2 conditions
                            beta22 = b*e/float(a)
                            #gamma2 conditions
                            g21 = c - e*e/float(c)
                            if e%c == 0 and beta22%c == 0 and g21%f == 0:
                                d = 0
                                #alpha1 condition
                                while (d<f):
                                    #beta1 condition
                                    beta11 = -a+b+d
                                    beta12 = -a+b-b*d/float(a)
                                    #gammma conditions
                                    g11 = -a+b+d-e*beta11/float(c)
                                    g12 = -a+b+d-d*d/float(a) -e*beta12/float(c)
                                    g22 = c-d*e/float(a) + e*beta22/float(c)
                                    if beta11%c == 0 and beta12%c==0 and g11%f==0 and g12%f==0 and g22%f==0:
                                        HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                        srHNFs.append(HNF)
                                    d +=a
                            e += a
                                
    return srHNFs
