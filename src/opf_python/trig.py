def trig_spHNFs(n):
    """Finds the symmetry preserving HNFs for the trigonal lattices
    with a determinant of n. Assuming A = [[1,2,2],[2,1,2],[4,3,3]].

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        spHNFs (list of lists): The symmetry preserving HNFs.
    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    spHNFs = []
    
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
                                    spHNFs.append(HNF)
                                
                                
                                d += a
                        e += a
    return spHNFs

def rhom_2_4(n):
    """Finds the symmetry preserving HNFs for the rhombohedral lattice
    that matches niggli conditions for niggli cells number 2 and 4
    with a determinant of n. Assuming A = [[-1, 0,-1],[0, -1.32288,
    -0.5],[-1.11652, -0.610985, 0.616515]] for number 2 and A= [[-1,
    0,-1],[0, -1.32288, 0.5],[-0.548584, 0.774292, 1.04858]].

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        spHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    spHNFs = []
    
    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]

        # conditions alpha 3 and betta3
        if f%a == 0 and f%c == 0:
            # bs = range(0,c)
            # es = range(0,f,a)
            # ds = range(0,f,a)
            
            if c<f:
                e = f-c # By inspecting results.
            else:
                e = 0
            bs = range(0,c)
            if e!=0:
                if a <c:
                    bs = [0,a]
                else:
                    bs = [0]
            else:
                if a==c:
                    bs = [0]
                else:
                    if a<c:
                        bs = [a]
                    else:
                        bs = [0]
            if c==1:
                ds = [e]
            elif e==0:
                if a==c:
                    ds = [0]
                else: 
                    ds = [a]
            else:
                ds = range(0,f,a)

            # for e in es:
            g12 = -c+e*e/float(c)
            #beta2 condition
            if e%c==0 and g12%f==0 and e%a==0:
                for b in bs:
                    b32 = b*f/float(a)
                    b22 = b*e/float(a)
                    if b32%c==0 and b22%c==0:
                        for d in ds:
                            b11 = b-d
                            b12 = -a+b*d/float(a)
                            g11 = -b+d-b11*e/float(c)
                            g12 = -b+d*d/float(a)-b12*e/float(c)
                            g22 = -c+d*e/float(a) -b*e*e/float(a*c)
                            if b11%c==0 and b12%c==0 and g11%f==0 and g12%f==0 and g22%f==0:
                                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                spHNFs.append(HNF)
                                
    return spHNFs

def rhom_9(n):
    """Finds the symmetry preserving HNFs for the rhombohedral lattice
    that matches niggli conditions for niggli cells number 2 and 4 with a
    determinant of n. Assuming A = [[0,-1,1],[1,-1,0],[2,1,2]].

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        spHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    spHNFs = []
    
    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]

        e = 0#By inspection of results.

        if f%a==0 and f%c==0 and c%a==0:
            for b in range(0,c,a):
                b31 = b*f/float(a)
                if b31%c==0:
                    if c ==1: # By inspection of results.
                        ds = [0]
                    else:
                        ds = range(0,f,a*f/c)
                    for d in ds:
                        a11 = -b-d
                        b11 = 2*b -b*a11/float(a)
                        b12 = a+2*b+d
                        g11 = 2*d -d*a11/float(a)
                        g21 = c*d/float(a)
                        if a11%a==0 and b11%c==0 and b12%c==0 and g11%f==0 and g21%f==0:
                            HNF = [[a,0,0],[b,c,0],[d,e,f]]
                            spHNFs.append(HNF)
                                
    return spHNFs                                    


def rhom_24(n):
    """Finds the symmetry preserving HNFs for the rhombohedral lattice
    that matches niggli conditions for niggli cells number 2 and 4
    with a determinant of n. Assuming A =
    [[-1,0,-1],[1.51184,0,-0.845178],[-0.255922,-1.44338,0.92259].

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        spHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    spHNFs = []
    
    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]
        b = 0
        
        if f%a==0 and f%c==0:
            # if c==1:
            #     es = [f-1]
            #     if f>3:
            #         ds = [3]
            #     else:
            #         ds = range(f)
            if c==f:
                es = [0]
                ds = [0]
            else: 
                es = [f-c]#range(0,f,a)
                if c==1 and f>3:
                    ds = [3]
                else:
                    ds = range(0,es[0]+1,c)
            for e in es:
                if e%c==0 and e%a==0:
                    for d in ds:
                        if d%c==0 and d%a==0:
                            # print("a",a,"b",b,"c",c,"d",d,"e",e,"f",f)
                            g11 = 2*d-d*d/float(a) -d*e/float(c)
                            g21 = 2*e-d*e/float(a)-e*e/float(c)
                            g22 = -c+e-d*e/float(a)-e*e/float(c)
                            # print("b11",b11,"b12",b12,"g11",g11,"g12",g12,"g21",g21,"g22",g22)
                            if g11%f==0 and g21%f==0 and g22%f==0:
                                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                spHNFs.append(HNF)
                                
    return spHNFs                                    
                                    
