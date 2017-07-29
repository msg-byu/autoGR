def body_ortho_srHNFs(n):
    """Finds the symmetry preserving HNFs for the body centered
    orthorhombic lattices with a determinant of n. Assuming A =
    [[0.5,1,1.5],[0,2,0],[0,0,3]] (3rd basis choince in
    notes/body_centered_ortho).

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        srHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    srHNFs = []
    
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
                                srHNFs.append(HNF)
                                
    return srHNFs
