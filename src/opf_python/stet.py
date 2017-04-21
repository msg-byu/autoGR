def stet_srHNFs(n):
    """Finds the symmetry preserving HNFs for the simple tetragonal lattices
    with a determinant of n. Assuming A = [[1,0,0],[0,1,0],[0,0,2]].

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
            #use beta1 condition (c divisible by 2) and gamma2
            #condition (f divisible by 2) to find b and e vals.
            if c%2==0:
                bs = [0,c/2]
            else:
                bs = [0]
            if f%2 ==0:
                es = [0,f/2]
            else:
                es = [0]
            #loop over values for e and b.
            for b in bs:
                #alpha1 condition
                if b%a==0:
                    beta13 = -a+b*b/float(a)
                    if beta13%c==0:
                        for e in es:
                            g12=2*b*e/float(c)
                            if g12%f==0:
                                for d in range(f):
                                    g13=-e*beta13/float(c)+d*(b/float(a)-1)
                                    g23=-e*(b/float(a)+1)+d*c/float(a)
                                    if g13%f==0 and g23%f==0:
                                        HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                        srHNFs.append(HNF)
                                
                                
    return srHNFs
