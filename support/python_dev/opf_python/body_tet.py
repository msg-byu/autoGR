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

def body_tet_6_7_15_18(n):
    """Finds the symmetry preserving HNFs for the body centered tetragonal
    lattices with a determinant of n. Assuming :
    A  = [[1.80278,-1.47253,0.762655],[2.80278,0.13535,-0.791285],
    [0.80278,-0.47253,2.762655]] for 6,
    A = [[1.95095, 1.19163, 0.879663],[0.0, 2.60788, 0.44606],
    [0.95095, -0.41625, 2.433603]] for 7,
    A = [[-1.0,-1.0,2.0],[0.0,-2.0,0.0],[-2.0,0.0,0.0]] for 15,
    A = [[-2.0,-1.0,1.0],[-3.0,1.0,0.0],[-1.0,-3.0,0.0]] for 18.


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

        if f%c == 0:
            if f%2==0:
                es = [0,f//2.]
            else:
                es = [0]
            for e in es:
                if e%c ==0:
                    g21 = -c +(e*e/float(c))
                    if g21%f == 0:
                        for d in range(f):
                            g13 = a+2*d
                            if g13%f == 0:
                                for b in range(c):
                                    b12 = b-d
                                    if b12%c == 0:
                                        g12 = -b+d-(b12*e/float(c))
                                        if g12%f == 0:
                                            HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                            spHNFs.append(HNF)

    return spHNFs

