def face_ortho_srHNFs(n):
    """Finds the symmetry preserving HNFs for the face centered
    orthorhombic lattices with a determinant of n. Assuming A =
    [[0,1,1.5],[0.5,0,1.5],[0,0,3]] (10th basis choince in
    notes/body_centered_ortho).

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

        #beta1 candition
        if c%2==0:
            bs = [0,c/2]
        else:
            bs = [0]

        for b in bs:
            for e in range(f):
                #gamma2 condition and gamma1 condition
                gamma23 = c+2*e
                gamma12 = b+2*b*e/float(c)
                if gamma23%f==0 and gamma12%f==0:
                    for d in range(f):
                        #gamma1 condition
                        gamma13 = a+b+2*d
                        if gamma13%f==0:
                            HNF = [[a,0,0],[b,c,0],[d,e,f]]
                            srHNFs.append(HNF)
                                
    return srHNFs
