def body_tet_15(n):
    """Finds the symmetry preserving HNFs for the body centered tetragonal
    lattices with a determinant of n. Assuming A =
    [[-1,1,2],[1,-1,2],[1,1,-2]] (first basis choince in notes/body_centered_tet.nb).

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        spHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    spHNFs = []
    primes = 0
    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]

        #alpha3 condition
        if a==c and a==f:
            bs = [0]
            ds = [0]
            es = [0]
        elif f==(a*c*f):
            bs = [0]
            if f%2==0:
                ds = [1,f//2 +1]
                es = [1,f//2 +1]
            else:
                ds = [1]
                es = [1]
        else:
            ds = range(0,f,smallest_prime(c))
            es = range(0,f,smallest_prime(c))
            if a==1:
                bs = [1,c//2+1]
            else:
                bs = range(0,c,a)

        #alpha3 condition
        if f%a==0:
            #beta3 condition
            if f%c ==0:
                for b in bs:
                    #beta3 condition
                    if (b*f)%(a*c) == 0 and b<c:
                        for e in es:
                            #beta2 conditions
                            beta22 = b*e/float(a)
                            #gamma2 conditions
                            g21 = c - e*e/float(c)
                            if e%c == 0 and beta22%c == 0 and g21%f == 0 and e<f and e%a==0:
                                for d in ds:
                                    #beta1 condition
                                    beta11 = -a+b+d
                                    beta12 = -a+b-b*d/float(a)
                                    #gammma conditions
                                    g11 = -a+b+d-e*beta11/float(c)
                                    g12 = -a+b+d-d*d/float(a) -e*beta12/float(c)
                                    g22 = c-d*e/float(a) + e*beta22/float(c)
                                    if beta11%c == 0 and beta12%c==0 and g11%f==0 and g12%f==0 and g22%f==0 and d<f and d%a==0:
                                        HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                        spHNFs.append(HNF)

    return spHNFs

def smallest_prime(a):
    """Finds the smallest prime factor of a.

    Args:
        a (int): The number we want the prime factor of.

    Return:
        pf (int): The smallest prime factor.
    """

    if a <= 2:
        return a
    else:
        for p in range(2,a+1):
            if a%p==0:
                return p
                                    
def body_tet_7(n):
    """Finds the symmetry preserving HNFs for the body centered tetragonal
    lattices with niggli setting 7 a determinant of n. Assuming A =
    [[-1.95095 , 1.41625 , -0.433603], [ 1.  , -1.  , -2.  ], [
    1.95095 , 1.19163 , 0.879663]].

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

        if f%a==0 and f%c==0 and a%c==0:
            if a==1 and c==1:
                es = [f-1]
                bs = [0]
            elif a==c and c==f:
                es = [0]
                bs = [0]
            else:
                es = [0]+list(range(c,f,a))
                bs = range(c)
            for b in bs:
                b11 = -a+2*b
                b33 = f-b*f/float(a)
                if b11%c==0 and b33%c==0:
                    for e in es:
                        a23 = -c+e
                        b23 = e-b*a23/float(a)
                        if a23%a==0 and b23%c==0:
                            if b==0:
                                if f%2==0 and (f//2+a)<f:
                                    ds = [0,a,f//2+a]
                                elif a<f:
                                    ds = [0,a]
                                else:
                                    ds = range(f)
                            else:
                                if e == f-c:
                                    ds = [f//2+b]
                                else:
                                    ds = []
                                    growing = True
                                    count = 0
                                    while growing:
                                        temp = (count*f+e*b11/float(c)+a)//2
                                        if temp < f and temp > 0:
                                            ds.append(temp)
                                        if temp >f:
                                            growing = False
                                        count += 1
                                
                            for d in ds:
                                a13 = -b+d
                                b13 = -a+d-b*a13/float(a)
                                g11 = -a+2*d -e*b11/float(c)
                                g13 = d-d*a13/float(a) - e*b13/float(c)
                                g23 = e-d*a23/float(a)-e*b23/float(c)
                                if a13%a==0 and b13%c==0 and g11%f==0 and g13%f==0 and g23%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    spHNFs.append(HNF)

    return spHNFs
                                    
def body_tet_6(n):
    """Finds the symmetry preserving HNFs for the body centered tetragonal
    lattices with niggli setting 6 a determinant of n. Assuming A =
    [[-1.  , 1.  , 2.  ], [ 1.  , 1.60788 , -1.55394 ], [ 1.80278 ,
    -1.47253 , 0.762655]].

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

        if f%c==0 and f%a==0:
            if a==c and f==a:
                es = [0]
                ds = [0]
                bs = [0]
            else:
                es = [0]+list(range(c,f,f//2))
                if c%2==0:
                    ds = range(0,f,c//2)
                else:
                    ds = range(0,f,c)
                bs = range(c)
            for b in bs:
                if (b*f)%(a*c)==0:
                    for e in es:
                        g21 = c-e*e/float(c)
                        b22 = b*e/float(a)
                        if e%a==0 and g21%f==0 and b22%c==0:
                            for d in ds:
                                b11 = -a+b+d
                                b12 = -a+b-b*d/float(a)
                                g11 = b11 - b11*e/float(c)
                                g12 = b11 -d*d/float(a) -b12*e/float(c)
                                g22 = c -d*e/float(a) +b22*e/float(c)
                                if b11%c==0 and b12%c==0 and g11%f==0 and g12%f==0 and g22%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    spHNFs.append(HNF)

    return spHNFs

def body_tet_18(n):
    """Finds the symmetry preserving HNFs for the body centered tetragonal
    lattices with niggli setting 18 a determinant of n. Assuming A =
    [[ 0, 0, 2], [ 1, -2, 1], [-2, -1, 1]].

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

        if c%a==0 and f%a==0 and f%c==0:
            if c==f:
                es = [0]
                if a==c:
                    bs = [0]
                    ds = [0]
                else:
                    ds = range(0,f,smallest_prime(c))
                    bs = range(0,c,smallest_prime(c))
            elif a*c*f==f:
                bs = [0]
                if f%2==0:
                    es = [f//2-1,f-1]
                else:
                    es = [f-1]
                if f>=2:
                    ds = [f-2]
                else: # pragma: no cover
                    ds = [0]
            else:
                bs = range(0,c,smallest_prime(c))
                if f%2==0:
                    es = [f//2-c,f-c]
                else:
                    es = [f-c]
                if c%2==0:
                    ds = range(0,f,c//2)
                else:
                    ds = range(0,f,c)
                
            for b in bs:
                b31 = f+b*f/float(a)
                b33 = b*f/float(a)
                if b31%c==0 and b33%c==0:
                    for e in es:
                        a21 = -c-e
                        b23 = -b*a21/float(a)
                        b21 = b23+e
                        b22 = b*c/float(a) -e
                        if a21%a==0 and b23%c==0 and b21%c==0 and b22%c==0 and b23%c==0:
                            for d in ds:
                                a11 = -b-d
                                b11 = b-b*a11/float(a)+d
                                b12 = b+b*b/float(a)-d
                                b13 = 2*b-b*a11/float(a)
                                g11 = b+d-a11*d/float(a)-b11*e/float(c)
                                g12 = b+d+b*d/float(a)-b12*e/float(c)
                                g13 = 2*d-a11*d/float(a)-b13*e/float(c)
                                g21 = c-d*a21/float(a)-e*b21/float(c)
                                g22 = c+c*d/float(a)-b22*e/float(c)
                                g23 = -d*a21/float(a)-b23*e/float(c)
                                if a11%a==0 and b11%c==0 and b12%c==0 and b13%c==0 and g11%f==0 and g12%f==0 and g13%f==0 and g21%f==0 and g22%f==0 and g23%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    spHNFs.append(HNF)

    return spHNFs
                                
