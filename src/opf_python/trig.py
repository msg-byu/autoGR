def trig_srHNFs(n):
    """Finds the symmetry preserving HNFs for the trigonal lattices
    with a determinant of n. Assuming A = [[1,2,2],[2,1,2],[4,3,3]].

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
            #beta1 condition
            if c%2==0:
                bs = [0,c/2]
            else:
                bs = [0]
                    
            for b in bs:
                #beta3 condition
                beta13 = f+b*f/float(a)
                if beta13%c==0:                    
                    #find e values alpha2 condition
                    e = 0
                    while (e<f):
                        beta22 = e+b*e/float(a)
                        g21= c+2*e
                        g11 = b+2*b*e/float(c)
                        if beta22%c==0 and g21%f==0 and g11%f==0:
                            #find d from alpha1 condition
                            d = 0
                            while (d<f):
                                #beta2 condition
                                beta12 = -a+b+d+d*b/float(a)
                                g12 = -b-d+d*d/float(a)-e*beta12/float(c)
                                g22 = -c-2*e+d*e/float(a)-e*beta22/float(c)
                                if beta12%c==0 and g12%f==0 and g22%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    srHNFs.append(HNF)
                                
                                
                                d += a
                        e += a
    return srHNFs
