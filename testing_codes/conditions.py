import numpy as np

def get_HNF_diagonals(n):
    """Finds the diagonals of the HNF that reach the target n value.
    
    Args:
        n (int): The target determinant for the HNF.
        
    Retruns:
        diags (list of lists): The allowed values of the determinant.
    """
    
    diags = []
    for i in range(1,n+1):
        if not n%i == 0:
            continue
        else:
            q = n/i
            for j in range(1,q+1):
                if not q%j == 0:
                    continue
                else:
                    diags.append([i,j,q/j])
                    
    return diags

def forms_group(gens,pg):
    """Tests if the given generators forms a group.
    
    Args:
        gens (list of list): The generators to check.
        pg (list of list): The group the generators form.
        
    Returns:
        corret_gens (bool): True if the generators form the group.
    """
    
    correct_gens = False
    
    group = []
    for i in gens:
        for j in gens:
            test = np.matmul(i,j)
            in_group = False
            for k in group:
                if np.allclose(test,k):
                    in_group = True

            if not in_group:
                group.append(test)
    
    growing = True
    while growing:
        nfound = 0
        for i in gens:
            for j in group:
                test = np.matmul(i,j)
                in_group = False
                for k in group:
                    if np.allclose(test,k):
                        in_group = True
                if not in_group:
                    group.append(test)
                    nfound += 1
                    
        if nfound == 0:
            growing = False

    if not len(pg) == len(group):
        correct_gens = False
        
    else:
        for i in pg:
            in_group = False
            for k in group:
                if np.allclose(i,k):
                    correct_gens = True
                    break
            if correct_gens == False:
                break

    return correct_gens

def find_gens_of_pg(pg):
    """This subroutine finds the generators of the point group.
    
    Args:
        pg (list of list): A list of the matrix form of the point group.
        
    Returns:
        gens (list of list): Those operations that will generate the 
            remainder of the group.
    """
    
    from itertools import combinations
    
    n_gens = 1
    found_gens = False
    while not found_gens:
        possible_gens = list(combinations(range(len(pg)),r=n_gens))
        for test in possible_gens:
            test_gens = []
            for i in test:
                test_gens.append(pg[i])            
            if forms_group(test_gens,pg):
                gens = test_gens
                found_gens = True
                break
        n_gens += 1

    return gens


def div_HNF(lat,n):
    """Finds the HNFs that preserve the symmetry of the lattice.
    
    Args:
        lat (numpy.ndarray): The vectors (as rows) of the parent lattice.
        n (int): The volume factor for the supercell.
        
    Returns:
        HNFs (list of lists): The HNFs the preserve the symmetry.
    """
    
    from phenum.symmetry import _get_lattice_pointGroup
    
    diags = get_HNF_diagonals(n)
    
    pg = _get_lattice_pointGroup(lat)
    gens = find_gens_of_pg(pg)

    # transpose the lattice so that it has the right form for the rest of the 
    # operations.
    lat = np.transpose(lat)

    lat_gens = []
    for g in gens:
        temp = np.matmul(np.linalg.inv(lat),np.matmul(g,lat))
        lat_gens.append(np.transpose(temp))
    
    x11 = []
    x12 = []
    x13 = []
    x21 = []
    x22 = []
    x23 = []
    x31 = []
    x32 = []
    x33 = []
    
    for g in lat_gens:
        # print("g",g)
        x11.append(g[0][0])
        x12.append(g[0][1])
        x13.append(g[0][2])
        x21.append(g[1][0])
        x22.append(g[1][1])
        x23.append(g[1][2])
        x31.append(g[2][0])
        x32.append(g[2][1])
        x33.append(g[2][2])
    
    x11 = np.array(x11)
    x12 = np.array(x12)
    x13 = np.array(x13)    
    x21 = np.array(x21)
    x22 = np.array(x22)
    x23 = np.array(x23) 
    x31 = np.array(x31)
    x32 = np.array(x32)
    x33 = np.array(x33)    
    
    count = 0
    HNFs = []

    for diag in diags:
        print("diag",diag)
        a = diag[0]
        c = diag[1]
        f = diag[2]
    
        # a divides tests
        if np.allclose((x13*f)%a,0):
            d = None
            e = None
            b = None
            if np.allclose(x13,0) and not np.allclose(x12,0):
                # d and e are unknown and b=0.
                if not np.allclose((x12*c)%a,0):
                    # print("c cond",(x12*c)%a)
                    continue
                b = 0
                al1 = b*x12/a
                al2 = c*x12/a
                al3 = f*x13/a
                tHNFs = cdivs(a,b,c,d,e,f,al1,al2,al3,x11,x21,x22,x23,x31,x32,x33)
                for t in tHNFs:
                    HNFs.append(t)
                    count += 1
            elif np.allclose(x12,0) and not np.allclose(x13,0):
                # b is unkown but d and e can have same values.
                vals = []
                N = 0
                xt = x13[np.nonzero(x13)]
                val = np.unique(N*a/xt)
                while any(abs(val) < f):
                    for v in val:
                        if v < f:
                            vals.append(v)
                    N += 1
                    val = np.unique(N*a/xt)
                for d in vals:
                    for e in vals:
                        al1 = d*x13/a
                        al2 = e*x13/a
                        al3 = f*x13/a
                        tHNFs = cdivs(a,b,c,d,e,f,al1,al2,al3,x11,x21,x22,x23,x31,x32,x33)
                        for t in tHNFs:
                            HNFs.append(t)
                            count += 1
            else:
                for e in range(f):
                    if np.allclose((c*x12 +e*x13)%a,0):
                        for b in range(c):
                            for d in range(f):
                                if np.allclose((b*x12+d*x13)%a,0):
                                    al1 = (b*x12+d*x13)/a
                                    al2 = (c*x12+e*x13)/a
                                    al3 = f*x13/a
                                    tHNFs = cdivs(a,b,c,d,e,f,al1,al2,al3,x11,x21,x22,x23,x31,x32,x33)
                                    for t in tHNFs:
                                        HNFs.append(t)
                                        count += 1
                                else:
                                    continue
                    else:
                        continue
                                
        else:
            # print("f cond")
            continue
    
    return HNFs

def fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33):
    """Finds the f divides conditions for the symmetry preserving HNFs.
    
    Args:
        a (int): a from the HNF.
        b (int): b from the HNF.
        c (int): c from the HNF.
        d (int): d from the HNF.
        e (int): e from the HNF.
        f (int): f from the HNF.
        al1 (numpy.array): array of alpha1 values from write up.
        al2 (numpy.array): array of alpha2 values from write up.
        be1 (numpy.array): array of beta1 values from write up.
        be2 (numpy.array): array of beta2 values from write up.
        x11 (numpy.array): array of pg values for x(1,1) spot.
        x22 (numpy.array): array of pg values for x(2,2) spot.
        x31 (numpy.array): array of pg values for x(3,1) spot.
        x32 (numpy.array): array of pg values for x(3,2) spot.
        x33 (numpy.array): array of pg values for x(3,3) spot.
        
    Returns:
        HNFs (list of lists): The symmetry preserving HNFs.
    """
    # print("***************enter fdivs")
    # print("b: ",b," d: ",d," e: ",e)
    HNFs = []
    
    if b == None and d == None and e == None:
        xvar1 = (x33-x22-be2)
        xvar2 = (x33-x11-al1)
        for b in range(c):
            for e in range(f):
                if not np.allclose(xvar2,0):
                    N = min(np.round((a*x31+b*x32-be1*e)/f))
                    xt = xvar2[np.nonzero(xvar2)]
                    val = np.unique(np.reshape(np.outer(N*f-a*x31-b*x32+be1*e,1/xt),len(xt)*len(x32)))
                    while any(abs(val)<f):
                        for v in val:
                            if v < f and v >= 0 and np.allclose(v%1,0):
                                d = v
                                f1 = a*x31+b*x32+d*var2-be1*e
                                f2 = c*x32-d*al2+e*(x33-x33-be2)
                                if np.allclose(f1%f,0) and np.allclose(f2%f,0):
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    HNFs.append(HNF)
                        N += 1
                        val = np.unique(np.reshape(np.outer(N*f-a*x31-b*x32+be1*e,1/xt),len(xt)*len(x32)))
                elif not np.allclose(al2,0):
                    N = max(np.round((c*x32+e*var1)/f))
                    at = al2[np.nonzero(al2)]
                    val = np.unique(np.reshape(np.outer(-N*f+c*x32+e*var1,1/at),len(x32)*len(at)))
                    while any(abs(val)<f):
                        for v in val:
                            if v < f and v >= 0 and np.allclose(v%1,0):
                                d = v
                                f1 = a*x31+b*x32+d*var2-be1*e
                                f2 = c*x32-d*al2+e*(x33-x33-be2)
                                if np.allclose(f1%f,0) and np.allclose(f2%f,0):
                                    HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                    HNFs.append(HNF)
                        N -= 1
                        val = np.unique(np.reshape(np.outer(-N*f+c*x32+e*var1,1/at),len(x32)*len(at)))
                else:
                    for d in range(f):
                        f1 = a*x31+b*x32+d*var2-be1*e
                        f2 = c*x32-d*al2+e*(x33-x33-be2)
                        if np.allclose(f1%f,0) and np.allclose(f2%f,0):
                            HNF = [[a,0,0],[b,c,0],[d,e,f]]
                            HNFs.append(HNF)
            
    elif b == None:
        f2 = c*x32-d*al2+e*(x33-x22-be2)
        if np.allclose(f2%f,0):
            if not np.allclose(x32,0):
                N = min(np.round(a*x31+d*(x33-x11-al1)/f))
                xt = x32[np.nonzero(x32)]
                val = np.unique(np.reshape(np.outer(N*f-a*x31-d*(x33-x11-al1),1/xt),len(x33)*len(xt)))
                while any(abs(val)<c):
                    for v in val:
                        if v<c and v>=0 and np.allclose(v%1,0):
                            b = v
                            f1 = a*x32 + b*x32 + e*be1 +d*(x33-x11-al1)
                            if np.allclose(f1%f,0):
                                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                HNFs.append(HNF)
                    N += 1
                    val = np.unique(np.reshape(np.outer(N*f-a*x31-d*(x33-x11-al1),1/xt),len(x33)*len(xt)))
            else:
                for b in range(c):
                    f1 = a*x32 + b*x32 + e*be1 +d*(x33-x11-al1)
                    if np.allclose(f1%f,0):
                        HNF = [[a,0,0],[b,c,0],[d,e,f]]
                        HNFs.append(HNF)
        
    elif d==None and e == None:
        for e in range(f):
            if not np.allclose(xvar2,0):
                N = min(np.round((a*x31+b*x32-be1*e)/f))
                xt = xvar2[np.nonzero(xvar2)]
                val = np.unique(np.reshape(np.outer(N*f-a*x31-b*x32+be1*e,1/xt),len(xt)*len(x32)))
                while any(abs(val)<f):
                    for v in val:
                        if v < f and v >= 0 and np.allclose(v%1,0):
                            d = v
                            f1 = a*x31+b*x32+d*var2-be1*e
                            f2 = c*x32-d*al2+e*(x33-x33-be2)
                            if np.allclose(f1%f,0) and np.allclose(f2%f,0):
                                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                HNFs.append(HNF)
                    N += 1
                    val = np.unique(np.reshape(np.outer(N*f-a*x31-b*x32+be1*e,1/xt),len(xt)*len(x32)))
            elif not np.allclose(al2,0):
                N = max(np.round((c*x32+e*var1)/f))
                at = al2[np.nonzero(al2)]
                val = np.unique(np.reshape(np.outer(-N*f+c*x32+e*var1,1/at),len(x32)*len(at)))
                while any(abs(val)<f):
                    for v in val:
                        if v < f and v >= 0 and np.allclose(v%1,0):
                            d = v
                            f1 = a*x31+b*x32+d*var2-be1*e
                            f2 = c*x32-d*al2+e*(x33-x33-be2)
                            if np.allclose(f1%f,0) and np.allclose(f2%f,0):
                                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                                HNFs.append(HNF)    
                    N -= 1
                    val = np.unique(np.reshape(np.outer(-N*f+c*x32+e*var1,1/at),len(x32)*len(at)))
            else:
                for d in range(f):
                    f1 = a*x31+b*x32+d*var2-be1*e
                    f2 = c*x32-d*al2+e*(x33-x33-be2)
                    if np.allclose(f1%f,0) and np.allclose(f2%f,0):
                        HNF = [[a,0,0],[b,c,0],[d,e,f]]
                        HNFs.append(HNF)
    else:
        if e==None or d==None or b==None:
            print("*****************ERROR IN fdivs**************")
        else:
            f2 = c*x32-d*al2+e*(x33-x22-be2)
            f1 = a*x31 + b*x32 + e*be1 +d*(x33-x11-al1)
           
            # if np.allclose(e,1) and np.allclose(d,1):
            #     print("e: ",e," d: ",d," b: ",b)
            #     print("x31: ",x31)
            #     print("x32: ",x32)
            #     print("al2: ",al2)
            #     print("x33: ",x33)
            #     print("x22: ",x22)
            #     print("be2: ",be2)
            #     print("al1: ",al1)
            #     print("f1: ",f1)
            #     print("f2: ",f2)
            if np.allclose(f1%f,0) and np.allclose(f2%f,0):
                HNF = [[a,0,0],[b,c,0],[d,e,f]]
                HNFs.append(HNF)
    # print("***********exit fdivs**************")
    return HNFs

def cdivs(a,b,c,d,e,f,al1,al2,al3,x11,x21,x22,x23,x31,x32,x33):
    """Finds the c divides conditions for the symmetry preserving HNFs.
    
    Args:
        a (int): a from the HNF.
        b (int): b from the HNF.
        c (int): c from the HNF.
        d (int): d from the HNF.
        e (int): e from the HNF.
        f (int): f from the HNF.
        al1 (numpy.array): array of alpha1 values from write up.
        al2 (numpy.array): array of alpha2 values from write up.
        al3 (numpy.array): array of alpha3 values from write up.
        x11 (numpy.array): array of pg values for x(1,1) spot.
        x21 (numpy.array): array of pg values for x(2,1) spot.
        x22 (numpy.array): array of pg values for x(2,2) spot.
        x23 (numpy.array): array of pg values for x(2,3) spot.
        x31 (numpy.array): array of pg values for x(3,1) spot.
        x32 (numpy.array): array of pg values for x(3,2) spot.
        x33 (numpy.array): array of pg values for x(3,3) spot.
        
    Returns:
        HNFs (list of lists): The symmetry preserving HNFs.
    """
    
    HNFs = []
    if np.allclose(x23,0):
        if b == None:
            # find the b values, d and e still unkown
            if not np.allclose(al3, 0):
                N=0
                at = al3[np.nonzero(al3)]
                val = np.unique(N*c/at)
                while any(abs(val) <c):
                    for v in val:
                        if v < c and v >= 0 and np.allclose(v%1==0):
                            b = v
                            c1 = a*x21 + b*(x22-al1-x11)
                            c2 =(-b*al2)
                            if np.allclose(c1%c,0) and np.allclose(c2%c,0):
                                be1 = c1/c
                                be2 =c2/c
                                tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33) 
                                for t in tHNFs:
                                    HNFs.append(t)
                    N += 1
                    val = np.unique(N*c/at)
            elif not np.allclose(al2,0):
                N=0
                at = al2[np.nonzero(al2)]
                val = np.unique(N*c/at)
                while any(abs(val) <c):
                    for v in val:
                        if v < c and v>=0 and np.allclose(v%1,0):
                            b = v
                            c1 = a*x21 + b*(x22-al1-x11)
                            c3 =(-b*al3)
                            if np.allclose(c1%c,0) and np.allclose(c3%c,0):
                                be1 = c1/c
                                be2 =-b*al2/c
                                tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)                        
                                for t in tHNFs:
                                    HNFs.append(t)
                    N += 1
                    val = np.unique(N*c/at)
            else:
                if not np.allclose((x22-x11-al1),0):
                    N=0
                    xt = (x22-x11-al1)
                    xt = xt[np.nonzero(xt)]
                    val = np.unique(np.reshape(np.outer(N*c-a*x21,1/xt),len(x21)*len(xt)))
                    while any(abs(val) <c):
                        for v in val:
                            if v < c and v>=0 and np.allclose(v%1,0):
                                b = v
                                c2 = -b*al2
                                c3 =(-b*al3)
                                if np.allclose(c2%c,0) and np.allclose(c3%c,0):
                                    be1 = (a*x21+b*(x22-x11-al1))/c
                                    be2 =-b*al2/c
                                    tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)                        
                                    for t in HNFs:
                                        HNFs.append(t)
                        N += 1
                        xt = (x22-x11-al1)
                        xt = xt[np.nonzero(xt)]                        
                        val = np.unique(np.reshape(np.outer(N*c-a*x21,1/xt),len(x21)*len(xt)))
                    else:
                        c1 = a*x21
                        c2 = 0
                        c3 = 0
                        if np.allclose(c1%c,0) and np.allclose(c2%c,0) and np.allclose(c3%c,0):
                            tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)                        
                            for t in HNFs:
                                HNFs.append(t)
                
        else:
            c1 = a*x21 + b*(x22-al1-x11)
            c2 = (-b*al2)
            c3 = (-b*a13)
            if np.allclose(c1%c,0) and np.allclose(c2%c,0) and np.allclose(c3%c,0):
                tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                for t in HNFs:
                    HNFs.append(t)
            
    else:
        if np.allclose(al3,0):
            if np.allclose((f*x23)%c,0):
                if b == None and e == None and d == None:
                    if np.allclose(al3,0) and np.allclose(al2,0) and np.allclose(al3,0):
                        N = 0
                        xt = x23[np.nonzero(x23)]
                        val = np.unique(N*c/xt)
                        while any(abs(val)<f):
                            for v in val:
                                if v <f and v>=0 and np.allclose(v%1,0):
                                    e = v
                                    for b in range(c):
                                        N2 =0
                                        xt = x23[np.nonzero(x23)]
                                        val2 = np.unique(np.reshape(np.outer((N2*c-a*x21-b*(x22-x11)),1/xt),len(x22)*len(xt)))
                                        while any(abs(val2)<f):
                                            for v2 in val2:
                                                if v2 <f and v2>=0 and np.allclose(v2%1,0):
                                                    d = v2
                                                    be1 = (a*x21+b*(x22-x11)+d*x23)/c
                                                    be2 = e*x23/c
                                                    tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                                                    for t in tHNFs:
                                                        HNFs.appned(t)
                                            N2 += 1
                                            xt = x23[np.nonzero(x23)]
                                            val2 = np.unique(np.reshape(np.outer((N2*c-a*x21-b*(x22-x11)),1/xt),len(x22)*len(xt)))
                                    
                            N += 1
                            val = np.unique(N*c/xt)

                    elif not np.allclose(al3,0):
                        N = max(np.round(f*x23/c))
                        at = al3[np.nonzero(al3)]
                        val = np.unique(np.reshape(np.outer(-N*c+f*x23,1/at),len(x23)*len(al3)))
                        while any(abs(val) < c):
                            for v in val:
                                if v < c and v>=0 and np.allclose(v%1,0):
                                    b = v
                                    N2 = min(np.round(-b*al2/c))
                                    xt = x23[np.nonzero(x23)]
                                    val2 = np.unique(np.reshape(np.outer(N2*c+b*al2,1/xt),len(xt)*len(al2)))
                                    while any(abs(val2)<f):
                                        for v2 in val2:
                                            if v2 <f and v2>=0 and np.allclose(v2%1,0):
                                                e = v2
                                                N3 = min(np.round((a*x21+b*(x22-x11-al1))/c))
                                                xt = x23[np.nonzero(x23)]
                                                val3 = np.unique(np.reshape(np.outer(N3*c-a*x21-b*(x22-x11-al1),1/xt),len(xt)*len(x22)))
                                                while any(abs(val2)<f):
                                                    for v3 in val3:
                                                        if v3 <f and v3>=0 and np.allclose(v3%1,0):
                                                            d = v3
                                                            be1 = (a*x21+b*(x22-x11-al1)+d*x23)/c
                                                            be2 = (e*x32-b*al2)/c
                                                            tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                                                            for t in tHNFs:
                                                                HNFs.append(t)
                                                    N3 += 1
                                                    xt = x23[np.nonzero(x23)]
                                                    val3 = np.unique(np.reshape(np.outer(N3*c-a*x21-b*(x22-x11-al1),1/xt),len(xt)*len(x22)))
                                        N2 += 1
                                        xt = x23[np.nonzero(x23)]
                                        val2 = np.unique(np.reshape(np.outer(N2*c+b*al2,1/xt),len(x22)*len(xt)))
                            N -= 1
                            at = al3[np.nonzero(al3)]
                            val = np.unique(np.reshape(np.outer(-N*c+f*x23,1/at),len(x23)*len(at)))
                    
                    else:
                        for b in range(c):
                            N2 = min(np.round(-b*al2/c))
                            xt = x23[np.nonzero(x23)]
                            val2 = np.unique(np.reshape(np.outer(N2*c+b*al2,1/xt),len(xt)*len(al2)))
                            while any(abs(val2)<f):
                                for v2 in val2:
                                    if v2 <f and v2 >= 0 and np.allclose(v2%1,0):
                                        e = v2
                                        N3 = min(np.round((a*x21+b*(x22-x11-al1))/c))
                                        xt = x23[np.nonzero(x23)]
                                        val3 = np.unique(np.reshape(np.outer(N3*c-a*x21-b*(x22-x11-al1),1/xt),len(x22)*len(xt)))
                                        while any(abs(val2)<f):
                                            for v3 in val3:
                                                if v3 <f and v3 >= 0 and np.allclose(v3%1,0):
                                                    d = v3
                                                    be1 = (a*x21+b*(x22-x11-al1)+d*x23)/c
                                                    be2 = (e*x32-b*al2)/c
                                                    tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                                                    for t in tHNFs:
                                                        HNFs.append(t)
                                            N3 += 1
                                            xt = x23[np.nonzero(x23)]
                                            val3 = np.unique(np.reshape(np.outer(N3*c-a*x21-b*(x22-x11-al1),1/xt),len(xt)*len(x22)))
                                N2 += 1
                                xt = x23[np.nonzero(x23)]
                                val2 = np.unique(np.reshape(np.outer(N2*c+b*al2,1/xt),len(al2)*len(xt)))
                
                elif b == None:
                    if not np.allclose(al3,0):
                        N = max(np.round(f*x23/c))
                        at = al3[np.nonzero(al3)]
                        val = np.unique(np.reshape(np.outer(-N*c+f*x23,1/at),len(x23)*len(at)))
                        while any(abs(val) < c):
                            for v in val:
                                if v < c and v>= 0 and np.allclose(v%1,0):
                                    b = v
                                    c1 = a*x21+b*(x22-x11-al1)+d*x23
                                    c2 = -b*al2+e*x23
                                    if np.allclose(c1%c,0) and np.allclose(c2%c,0):
                                        be1 = c1/c
                                        be2 = c2/c
                                        tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                                        for t in tHNFs:
                                            HNFs.append(t)
                            N -= 1
                            at = al3[np.nonzero(al3)]
                            val = np.unique(np.reshape(np.outer(-N*c+f*x23,1/at),len(x23)*len(at)))
                        
                    elif not np.allclose(al2,0):
                        N = max(np.round(e*x23/c))
                        at = al2[np.nonzero(al2)]
                        val = np.unique(np.reshape(np.outer(-N*c+e*x23,1/at),len(x23)*len(at)))
                        while any(abs(val) < c):
                            for v in val:
                                if v < c and v>= 0 and np.allclose(v%1,0):
                                    b = v
                                    c1 = a*x21+b*(x22-x11-al1)+d*x23
                                    c2 = -b*al2+e*x23
                                    if np.allclose(c1%c,0) and np.allclose(c2%c,0):
                                        be1 = c1/c
                                        be2 = c2/c
                                        tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                                        for t in tHNFs:
                                            HNFs.append(t)
                            N -= 1
                            at = al2[np.nonzero(al2)]
                            val = np.unique(np.reshape(np.outer(-N*c+e*x23,1/at),len(x23)*len(at)))
                    
                    else:
                        if not np.allclose((x22-x11-al1),0):
                            N = min(np.round((a*x21-d*x23)/c))
                            xt = (x22-x11-al1)
                            xt = xt[np.nonzero(xt)]
                            val = np.unique(np.reshape(np.outer(N*c-a*x21sd*x23,1/xt),len(x23)*len(xt)))
                            while any(abs(val) < c):
                                for v in val:
                                    if v < c and v>=0 and np.allclose(v%1,0):
                                        b = v
                                        c1 = a*x21+b*(x22-x11-al1)+d*x23
                                        c2 = -b*al2+e*x23
                                        if np.allclose(c1%c,0) and np.allclose(c2%c,0):
                                            be1 = c1/c
                                            be2 = c2/c
                                            tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                                            for t in tHNFs:
                                                HNFs.append(t)
                                N += 1
                                xt = (x22-x11-al1)
                                xt = xt[np.nonzero(xt)]                               
                                val = np.unique(np.reshape(np.outer(N*c-a*x21sd*x23,1/xt),len(x23)*len(xt)))
                        else:
                            c1 = a*x21+d*x23
                            c2 = e*x23
                            c3 = f*x23
                            if np.allclose(c1%c,0) and np.allclose(c2%c,0) and np.allclose(c3%c,0):
                                tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33) 
                                for t in tHNFs:
                                    HNFs.append(t)
                    
                elif d == None and e == None:
                    N2 = min(np.round(-b*al2/c))
                    xt = x23[np.nonzero(x23)]
                    val2 = np.unique(np.reshape(np.outer(N2*c+b*al2,1/xt),len(xt)*len(al2)))
                    while any(abs(val2)<f):
                        for v2 in val2:
                            if v2 <f and v2>=0 and np.allclose(v2%1,0):
                                e = v2
                                N3 = min(np.round((a*x21+b*(x22-x11-al1))/c))
                                xt = x23[np.nonzero(x23)]
                                val3 = np.unique(np.reshape(np.outer(N3*c-a*x21-b*(x22-x11-al1),1/xt),len(x22)*len(xt)))
                                while any(abs(val3)<f):
                                    for v3 in val3:
                                        if v3 <f and v3>=0 and np.allclose(v3%1,0):
                                            d = v3
                                            be1 = (a*x21+b*(x22-x11-al1)+d*x23)/c
                                            be2 = (e*x32-b*al2)/c
                                            tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                                            for t in tHNFs:
                                                HNFs.append(t)
                                    N3 += 1
                                    xt = x23[np.nonzero(x23)]
                                    val3 = np.unique(np.reshape(np.outer(N3*c-a*x21-b*(x22-x11-al1),1/xt),len(x22)*len(xt)))
                        N2 += 1
                        xt = x23[np.nonzero(x23)]
                        val2 = np.unique(np.reshape(np.outer(N2*c+b*al2,1/xt),len(xt)*len(al2)))
                
                else:
                    c1 = a*x21+b*(x22-al1-x11)+d*x23
                    c2 = -b*al2+e*x23
                    c3 = -b*al3+f*x23
                    if np.allclose(c1%c,0) and np.allclose(c2%c,0) and np.allclose(c3%c,0):
                        be1 = c1/c
                        be2 = c2/c
                        tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                        for t in tHNFs:
                            HNFs.append(t)
            # else:
                # print("f: ",f)
                # print("c: ",c)
                # print("x32: ",x32)
                # print("failed f*x32/c")
        else:
            if b==None and d==None and e==None:
                N = max(np.round(f*x23/c))
                at = al3[np.nonzero(al3)]
                val = np.unique(np.reshape(np.outer(-N*c+f*x23,1/at),len(x23)*len(at)))
                while any(abs(val) < c):
                    for v in val:
                        if v < c and v>= 0 and np.allclose(v%1,0):
                            b = v
                            N2 = min(np.round(-b*al2/c))
                            xt = x23[np.nonzero(x23)]
                            val2 = np.unique(np.reshape(np.outer(N2*c+b*al2,1/xt),len(xt)*len(al2)))
                            while any(abs(val2)<f):
                                for v2 in val2:
                                    if v2 <f and v2>=0 and np.allclose(v2%1,0):
                                        e = v2
                                        N3 = min(np.round((a*x21+b*(x22-x11-al1))/c))
                                        xt = x23[np.nonzero(x23)]
                                        val3 = np.unique(np.reshape(np.outer(N3*c-a*x21-b*(x22-x11-al1),1/xt),len(x22)*len(xt)))
                                        while any(abs(val3)<f):
                                            for v3 in val3:
                                                if v3 <f and v3>=0 and np.allclose(v3%1,0):
                                                    d = v3
                                                    c1 = a*x21+b*(x22-x11-al1)+d*x23
                                                    c2 = -b*al2+e*x23
                                                    if np.allclose(c1%c,0) and np.allclose(c2%c,0):
                                                        be1 = c1/c
                                                        be2 = c2/c
                                                        tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                                                        for t in tHNFs:
                                                            HNFs.append(t)
                                            N3 += 1
                                            xt = x23[np.nonzero(x23)]
                                            val3 = np.unique(np.reshape(np.outer(N3*c-a*x21-b*(x22-x11-al1),1/xt),len(x22)*len(xt)))
                                N2 += 1
                                xt = x23[np.nonzero(x23)]
                                val2 = np.unique(np.reshape(np.outer(N2*c+b*al2,1/xt),len(xt)*len(al2)))
                    N -= 1
                    at = al3[np.nonzero(al3)]
                    val = np.unique(np.reshape(np.outer(-N*c+f*x23,1/at),len(x23)*len(at)))
                
            elif b==None:
                N = max(np.round(f*x23/c))
                at = al3[np.nonzero(al3)]
                val = np.unique(np.reshape(np.outer(-N*c+f*x23,1/at),len(x23)*len(at)))
                while any(abs(val) < c):
                    for v in val:
                        if v < c and v>= 0 and np.allclose(v%1,0):
                            b = v
                            c1 = a*x21+b*(x22-x11-al1)+d*x23
                            c2 = -b*al2+e*x23
                            if np.allclose(c1%c,0) and np.allclose(c2%c,0):
                                be1 = c1/c
                                be2 = c2/c
                                tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                                for t in tHNFs:
                                    HNFs.append(t)
                    N -= 1
                    at = al3[np.nonzero(al3)]
                    val = np.unique(np.reshape(np.outer(-N*c+f*x23,1/at),len(x23)*len(at)))

            elif d==None and e==None:
                N2 = min(np.round(-b*al2/c))
                xt = x23[np.nonzero(x23)]
                val2 = np.unique(np.reshape(np.outer(N2*c+b*al2,1/xt),len(xt)*len(al2)))
                while any(abs(val2)<f):
                    for v2 in val2:
                        if v2 <f and v2>=0 and np.allclose(v2%1,0):
                            e = v2
                            N3 = min(np.round((a*x21+b*(x22-x11-al1))/c))
                            xt = x23[np.nonzero(x23)]
                            val3 = np.unique(np.reshape(np.outer(N3*c-a*x21-b*(x22-x11-al1),1/xt),len(x22)*len(xt)))
                            while any(abs(val3)<f):
                                for v3 in val3:
                                    if v3 <f and v3>=0 and np.allclose(v3%1,0):
                                        d = v3
                                        c1 = a*x21+b*(x22-x11-al1)+d*x23
                                        c2 = -b*al2+e*x23
                                        if np.allclose(c1%c,0) and np.allclose(c2%c,0):
                                            be1 = c1/c
                                            be2 = c2/c
                                            tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                                            for t in tHNFs:
                                                HNFs.append(t)
                                    N3 += 1
                                    xt = x23[np.nonzero(x23)]
                                    val3 = np.unique(np.reshape(np.outer(N3*c-a*x21-b*(x22-x11-al1),1/xt),len(x22)*len(xt)))
                        N2 += 1
                        xt = x23[np.nonzero(x23)]
                        val2 = np.unique(np.reshape(np.outer(N2*c+b*al2,1/xt),len(xt)*len(al2)))

            else:
                be1 = c1/c
                be2 = c2/c
                tHNFs = fdivs(a,b,c,d,e,f,al1,al2,be1,be2,x11,x22,x31,x32,x33)
                for t in tHNFs:
                    HNFs.append(t)
                
    return HNFs
