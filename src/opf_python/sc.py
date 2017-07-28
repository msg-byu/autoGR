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


        if a==c and c==f:
            b = 0
            d = 0
            e = 0
            HNF = [[a,0,0],[b,c,0],[d,e,f]]
            srHNFs.append(HNF)
            break
        elif a==c and f/float(a)==2:
            b = 0
            d = a
            e = a
            HNF = [[a,0,0],[b,c,0],[d,e,f]]
            srHNFs.append(HNF)
            break
        elif c==f and f/float(a)==2:
            b = a
            d = a
            e = 0
            HNF = [[a,0,0],[b,c,0],[d,e,f]]
            srHNFs.append(HNF)
            break
        
    return srHNFs
