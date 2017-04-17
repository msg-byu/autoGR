def bcc_srHNFs(n):
    """Finds the symmetry preserving HNFs for the base centered cubic lattices
    with a determinant of n. Assuming A = [[-1,1,1],[1,-1,1],[1,1,-1]]


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

        # print("a: ",a," c: ",c," f: ",f)
        #alpha3 condition
        if f%a==0:
            # print("past f%a")
            #beta3 condition
            if f%c==0:
                # print("past f%c")
                #find b values beta2 condition
                b = 0
                while (b<c):
                    # print("b: ",b)
                    #find e values alpha2 condition
                    e = 0
                    while (e<f):
                        # print("e :",e)
                        #beta2 condition
                        g22 = -c +e*e/float(c)
                        if e%c==0 and g22%f ==0:
                            # print("past e%c")
                            d = 0
                            #alpha1 condition
                            while (d<f):
                                # print("d: ",d)
                                beta11 = -a+b-b*d/float(a)
                                beta12 = b-d
                                g11 = -a+b-d*(d/float(a)-1)-e*beta11/float(c)
                                g12 = -b+d-beta12*e/float(c)
                                g21 = c-e*e*b/float(a*c) -d*e/float(a)
                                # print("b11",beta11)
                                # print("b12",beta12)
                                # print("g11",g11)
                                # print("g12",g12)
                                # print("g21",g21)
                                if beta11%c==0 and beta12%c==0 and g11%f==0 and g12%f==0 and g21%f==0:
                                    # print("past all")
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    srHNFs.append(HNF)
                                
                                
                                d += a
                        e += a
                    b += a
    return srHNFs
