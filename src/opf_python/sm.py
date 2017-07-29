def sm_srHNFs(n):
    """Finds the symmetry preserving HNFs for the simple monoclinic
    lattices with a determinant of n. Assuming A =
    [[2,0,0],[0,2,0],[0.5,0,2]].

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

        #beta1 candition
        if c%2==0:
            bs = [0,c/2]
        else:
            bs = [0]
        #gamma2 condition
        if f%2==0:
            es = [0,f/2]
        else:
            es = [0]
            
        for b in bs:
            for e in es:
                #gamma1 condition and gamma1 condition
                gamma12 = 2*b*e/float(c)
                if gamma12%f==0:
                    for d in range(f):
                        HNF = [[a,0,0],[b,c,0],[d,e,f]]
                        srHNFs.append(HNF)
                                
    return srHNFs
