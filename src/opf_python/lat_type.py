def lat_type(lat):
    """Finds the lattice type for the provided lattice vectors.
    
    Args:
        lat (numpy.array): the lattice vectors as columns of the matrix.
        
    Returns:
        latType (str): the type of lattice, i.e., sc, bcc, fcc, st, so....
        latBasis (numpy.array): The canonical basis for this lattice.
    """
    
    from phenum.vector_utils import _minkowski_reduce_basis
    from phenum.symmetry import _get_lattice_pointGroup
    
    clat = np.array(_minkowski_reduce_basis(lat,1E-10))
    #clat = lat
    
    latType = None
    latBasis = None
    
    a1 = clat[0]
    a2 = clat[1]
    a3 = clat[2]
    
    a1n = np.sqrt(np.dot(a1,a1))
    a2n = np.sqrt(np.dot(a2,a2))
    a3n = np.sqrt(np.dot(a3,a3))
    a1ta2 = np.dot(a1,a2)/abs(float(a1n*a2n))
    a1ta3 = np.dot(a1,a3)/abs(float(a1n*a3n))
    a2ta3 = np.dot(a2,a3)/abs(float(a2n*a3n))
    
    tiny = 1E-6
    
    p_count = 0
    if abs(a1ta2)<tiny:
        p_count += 1
    if abs(a1ta3)<tiny:
        p_count += 1
    if abs(a2ta3)<tiny:
        p_count += 1
    
    ## find pg to determine crystial family
    ## once family is found go through possible lattice system
    fam = len(_get_lattice_pointGroup(clat))
    if fam == 48: #cubic
        if p_count==3:
            latType = 'sc'
            latBasis = np.array([[1,0,0],[0,1,0],[0,0,1]])
        elif abs(a1ta2 +1./3.) <tiny or abs(a1ta3+1./3.)<tiny or abs(a2ta3+1./3.)<tiny:
            latType = 'bcc'
            latBasis = np.array([[-1,1,1],[1,-1,1],[1,1,-1]])
        elif abs(a1ta2-1./2.)<tiny or abs(a1ta3-1./2.)<tiny or abs(a2ta3-1./2.)<tiny:
            latType = 'fcc'
            latBasis = np.array([[0,1,1],[1,0,1],[1,1,0]])
        else:
            print("Could not identify lattice (cubic).")
    elif fam == 24: #hex
        latType = 'hex'
        latBasis = np.array([[1,0,0],[0.5,-0.8660254037844386,0],[0,0,2]])
    elif fam == 12: #trig
        latType = 'trig'
        latBasis = np.array([[1,2,2],[2,1,2],[4,3,3]])
    elif fam == 16: #tet
        if p_count==3:
            latType = 'stet'
            latBasis = np.array([[1,0,0],[0,1,0],[0,0,2]])
        else:
            latTyp = 'btet'
            latBasis = np.array([[-1,1,2],[1,-1,2],[1,1,-2]])
    elif fam == 8: #ortho
        if p_count == 3:
            latType = 'so'
            latBasis = np.array([[1,0,0],[0,2,0],[0,0,3]])
        elif p_count == 2:
            if ((abs(abs(a1ta2)-0.5)<tiny and abs(a1n-a2n)<tiny) or 
            (abs(abs(a1ta3)-0.5)<tiny and abs(a1n-a3n)<tiny)or 
            (abs(abs(a2ta3)-0.5)<tiny and abs(a2n-a3n)<tiny)):
                latType = 'hex'
                latBasis = np.array([[1,0,0],[0.5,-0.8660254037844386,0],[0,0,2]])
            else:
                latType = 'co'
                latBasis = np.array([[0.5,1,0],[0.5,-1,0],[0,0,3]])
        elif (p_count == 1 ):
            if abs(a1ta2)<tiny:
                pa3ta1n = np.linalg.norm(np.dot(a3,a1)/np.dot(a1,a1))
                pa3ta2n = np.linalg.norm(np.dot(a3,a2)/np.dot(a2,a2))
                if (pa3ta1n-a1n/2.)<tiny and (pa3ta2n -a2n/2.)<tiny:
                    latType = 'bo'
                    latBasis = np.array([[0.5,1,1.5],[0,2,0],[0,0,3]])
                else:
                    latType = 'fo'
                    latBasis = np.array([[0,1,1.5],[0.5,0,1.5],[0,0,3]])
            elif abs(a1ta3)<tiny:
                pa2ta1n = np.linalg.norm(np.dot(a2,a1)/np.dot(a1,a1))
                pa2ta3n = np.linalg.norm(np.dot(a2,a3)/np.dot(a3,a3))
                if (pa2ta1n-a1n/2.)<tiny and (pa2ta3n -a3n/2.)<tiny:
                    latType = 'bo'
                    latBasis = np.array([[0.5,1,1.5],[0,2,0],[0,0,3]])
                else:
                    latType = 'fo'                
                    latBasis = np.array([[0,1,1.5],[0.5,0,1.5],[0,0,3]])
            else:
                pa1ta2n = np.linalg.norm(np.dot(a1,a2)/np.dot(a2,a2))
                pa1ta3n = np.linalg.norm(np.dot(a1,a3)/np.dot(a3,a3))
                if (pa1ta2n-a2n/2.)<tiny and (pa1ta3n -a3n/2.)<tiny:
                    latType = 'bo'
                    latBasis = np.array([[0.5,1,1.5],[0,2,0],[0,0,3]])
                else:
                    latType = 'fo'
                    latBasis = np.array([[0,1,1.5],[0.5,0,1.5],[0,0,3]])
        elif(abs(a1n-a2n)<tiny and abs(a2n-a3n)<tiny):
            latType = 'bo'
            latBasis = np.array([[0.5,1,1.5],[0,2,0],[0,0,3]])
        elif (abs(a1n-a2n)>tiny and abs(a1n-a3n)>tiny and abs(a2n-a3n)>tiny):
            latType = 'fo'
            latBasis = np.array([[0,1,1.5],[0.5,0,1.5],[0,0,3]])
        else:
            print("Could not identify lattice (ortho)")
    elif fam == 4: #mono
        if p_count==2:
            latType = 'sm'
            latBasis = np.array([[2,0,0],[0,2,0],[0.5,0,2]])
        else:
            latType = 'cm'
            latBasis = np.array([[1,1,0],[0,2,0],[0.5,0,2]])
    elif fam == 2: #tric
        latType = 'tric'
        latBasis = None
    else:
        print("Could not indentify lattice.")


    return latType, latBasis
