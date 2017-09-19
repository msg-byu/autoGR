def body_tet_srHNFs(n):
    """Finds the symmetry preserving HNFs for the body centered tetragonal
    lattices with a determinant of n. Assuming A =
    [[-1,1,2],[1,-1,2],[1,1,-2]] (first basis choince in notes/body_centered_tet.nb).

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        srHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    srHNFs = []
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
                                        srHNFs.append(HNF)

    return srHNFs

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
    
def body_tet_srHNFs_2(n):
    """Finds the symmetry preserving HNFs for the body centered tetragonal
    lattices with a determinant of n. Assuming A =
    [[2,0,0],[0,2,0],[1,1,2]] (second basis choince in notes/body_centered_tet.nb).

    Args:
        n (int): The determinant of the HNFs.

    Returns:
        srHNFs (list of lists): The symmetry preserving HNFs.

    """

    from opf_python.universal import get_HNF_diagonals

    diags = get_HNF_diagonals(n)

    srHNFs = []
    primes = 0
    for diag in diags:
        a = diag[0]
        c = diag[1]
        f = diag[2]

        #infered conditions
        if a==c and a==f:
            bs = [0]
            ds = [0]
            es = [0]
        elif f==(a*c*f):
            bs = [0]
            # if (f//2)%2==0:
            if f%2 == 0:
                ds = [0,f//2]
                es = [0,f//2]
            else:
                ds = [0]
                es = [0]
        else:
            ds = list(range(0,f,a))
            bs = list(range(0,c))
            es = list(range(0,f,a))

        #alpha3 condition
        if f%a==0:
            for e in es:
                #alpha2  conditions
                if (e%a==0) and ((c+e)%a==0):
                    a23 = (-c-e)
                    for b in bs:
                        b23 = -(b*a23/a) -e
                        b33 = -f+(b*f)/a
                        #beta3 and beta2
                        if ((b*f)%(a*c) == 0) and ((b*e)%(a*c)==0) and (b23%c ==0) and (b33%c ==0):

                            for d in ds:
                                a13 = -b-d
                                b12 = 2*b+(b*d)/a
                                b13 = -a -(b*a13)/a -d
                                g12 = 2*d + d*d/a - (b12*e/c)
                                g13 = d-a13*d/a - b13*e/c
                                g22 = d*e/a - b*e*e/(a*c)
                                g23 = -d*a23/a + e - b23*e/c
                                if a13%a == 0 and b12%c == 0 and b13%c == 0 and g12%f==0 and g13%f==0 and g22%f==0 and g23%f==0:
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    srHNFs.append(HNF)

    return srHNFs
