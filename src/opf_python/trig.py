def trig_srHNFs(n):
    """Finds the symmetry preserving HNFs for the trigonal lattices
    with a determinant of n. Assuming A = [[1,2,2],[2,1,2],[2,2,1]].

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
                #find d values alpha1 condition
                d = 0
                while (d<f):
                    #find b from beta1 condition
                    N = -(d/c)
                    b = c*N+d
                    while (b<c):
                        if b >= 0:
                            beta12 = -a+b*d/float(a)
                            if beta12%c==0 and (b*f)%(a*c)==0:
                                #find e from alpha2 condition
                                e = 0
                                while (e<f):
                                    #beta2 condition
                                    g21 = -c+e*e/float(c)
                                    if e%c==0 and g21%f==0 and (b*e)%(a*c)==0:
                                        g11 = -b-e*(b-d)/float(c)+d
                                        g12 = -b-e*beta12/float(c)+d*d/float(a)
                                        g22 = -c+d*e/float(a) - b*e*e/float(a*c)
                                        if g11%f==0 and g12%f==0 and g22%f==0:
                                            HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                            srHNFs.append(HNF)
                                
                                
                                    e += a
                        N -= 1
                        b = -c*N+d
                    d += a
    return srHNFs
