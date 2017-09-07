def base_mono_srHNFs(n):
    """Finds the symmetry preserving HNFs for the base centered monoclinic
    lattices with a determinant of n. Assuming A =
    [[1,1,0],[0,2,0],[0.5,0,2]] (2nd basis choice in
    notes/base_centered_monoclinic.nb).

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

        #gamma2 candition
        if f%2==0:
            es = [0,f/2]
        else:
            es = [0]

        for b in range(c):
            #beta1 condition
            beta12 = a+2*b
            if beta12%c==0:
                for e in es:
                    #gamma1 condition
                    gamma12 = e*beta12/float(c)
                    if gamma12%f==0:
                        for d in range(f):
                            HNF = [[a,0,0],[b,c,0],[d,e,f]]
                            srHNFs.append(HNF)
                                
    return srHNFs
