def so_32(n):
    """Finds the symmetry preserving HNFs for the simple orthorhombic lattices
    with a determinant of n. Assuming A = [[1,0,0],[0,2,0],[0,0,3]].

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
        #gamma1 and gamma2 condition
        if f%2==0:
            ed_vals = [0,f/2]
        else:
            ed_vals = [0]

        for b in bs:
            for e in ed_vals:
                #gamma1 condition
                if (2*b*e)%(f*c)==0:
                    for d in ed_vals:
                        HNF = [[a,0,0],[b,c,0],[d,e,f]]
                        srHNFs.append(HNF)
                                
                                
    return srHNFs
