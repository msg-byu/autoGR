def base_ortho_srHNFs(n):
    """Finds the symmetry preserving HNFs for the base centered
    orthorhombic lattices with a determinant of n. Assuming A =
    [[0.5,1,0],[0.5,-1,0],[0,0,3]] (1st basis choince in
    notes/base_ortho).

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
                                srHNFs.append(HNF)
                b += a
                                
                                
    return srHNFs
