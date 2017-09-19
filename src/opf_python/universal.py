"""Subtroutines needed by all the systems."""

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
            q = n//i
            for j in range(1,q+1):
                if not q%j == 0:
                    continue
                else:
                    diags.append([i,j,q//j])
                    
    return diags

def find_srBs(A,kpt,exact=False):
    """Generates the k-point grid vectors for the target k-point density
    for the given lattice.

    Args:
        A (numpy.array): The lattice vectors as columns of a matrix.
        kpt (int): The target k-point density desired.
        exact (bool): True if only the target k-point density is to be used.

    Returns:
        Bs (numpy.array): list of symmetry preserving supercells.

    Raises:
        ValueError if the lattice type can't be identified.
    """

    from opf_python.lat_type import lat_type
    import numpy as np
    
    name, basis, order = lat_type(np.transpose(A))
    # name_w, basis = lat_type(np.transpose(A))
    # name, basis_w = lat_type(np.linalg.inv(A))

    # # We need to find the difference in volume between the canonical
    # # and the users basis.
    # rat = (np.linalg.det(A)/np.linalg.det(basis))**(1/3.)
    # basis = rat*basis
    
    # Find transformation matrix bteween basis.
    M = np.dot(basis,np.linalg.inv(A))

    if not exact:
        ns, mult = find_volumes(name,kpt)
    else:
        ns, mult = [kpt], 1.0

    srHNFs = []
    if name == "sc":
        from opf_python.sc import sc_srHNFs
        for n in ns:
            temp = sc_srHNFs(n)
            for t in temp:
                srHNFs.append(t)
        
    elif name == "bcc":
        from opf_python.bcc import bcc_srHNFs
        for n in ns:
            temp = bcc_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
            
    elif name == "fcc":
        from opf_python.fcc import fcc_srHNFs
        for n in ns:
            temp = fcc_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
        
    elif name == "hex":
        from opf_python.hx import hex_srHNFs
        for n in ns:
            temp = hex_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
        
    elif name == "trig":
        from opf_python.trig import trig_srHNFs
        for n in ns:
            temp = trig_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
        
    elif name == "stet":
        from opf_python.stet import stet_srHNFs
        for n in ns:
            temp = stet_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
       
    elif name == "btet":
        from opf_python.body_tet import body_tet_srHNFs
        for n in ns:
            temp = body_tet_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
        
    elif name == "so":
        from opf_python.so import so_srHNFs
        for n in ns:
            temp = so_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
        
    elif name == "co":
        from opf_python.base_ortho import base_ortho_srHNFs
        for n in ns:
            temp = base_ortho_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
        
    elif name == "bo":
        from opf_python.body_ortho import body_ortho_srHNFs
        for n in ns:
            temp = body_ortho_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
        
    elif name == "fo":
        from opf_python.face_ortho import face_ortho_srHNFs
        for n in ns:
            temp = face_ortho_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
        
    elif name == "sm":
        from opf_python.sm import sm_srHNFs
        for n in ns:
            temp = sm_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
        
    elif name == "cm":
        from opf_python.base_mono import base_mono_srHNFs
        for n in ns:
            temp = base_mono_srHNFs(n)
            for t in temp:
                srHNFs.append(t)        
        
    elif name == "tric":
        for n in ns:
            diags = get_HNF_diagonals(n)
            for diag in diags:
                a = diag[0]
                c = diag[1]
                f = diaf[2]
                for b in range(c):
                    for d in range(f):
                        for e in range(f):
                            HNF = np.array([[a,0,0],[b,c,0],[d,e,f]])
                            srHNFs.append(HNF*mult)
        
    else:
        raise ValueError("ERROR: unrecognized lattice type: ",name)

    Bs = []
    for H in srHNFs:
        B = np.dot(np.linalg.inv(M),np.dot(basis,H))
        Bs.append(B)

    return Bs, srHNFs
        
def find_volumes(lat_type,kpd):
    """Finds the allowed n's for a given lattice around the desired k-point density.
    
    Args:
        lat_type (str): The type of lattice.
        kpd (int): the target valume factor.

    Returns:
        ns (list of int): The target valume's to consider.
        mult (int): multiplication factor for triclinic.

    Raises:
        ValueError if the lattice type can't be identified.
    """

    from math import floor, ceil

    mult = None
    ns = []
    
    nt = float(kpd)
    if lat_type == "sc":
        nb = int(floor(nt**(1/3))-1)
        nmin = kpd-1000
        nmax = kpd+1000
        while nb**3<nmax:
            nc = nb**3
            for m in [1,2,4]:
                if m*nc>nmin and m*nc<nmax:
                    ns.append(m*nc)
            nb += 1   
        
    elif lat_type == "bcc":
        nb = int(floor(nt**(1/3))-1)
        nmin = kpd-1000
        nmax = kpd+1000
        while nb**3<nmax:
            nc = nb**3
            for m in [1,2,4]:
                if m*nc>nmin and m*nc<nmax:
                    ns.append(m*nc)
            nb += 1   
            
    elif lat_type == "fcc":
        nb = int(floor(nt**(1/3))-1)
        nmin = kpd-1000
        nmax = kpd+1000
        while nb**3<nmax:
            nc = nb**3
            for m in [1,4,16]:
                if m*nc>nmin and m*nc<nmax:
                    ns.append(m*nc)
            nb += 1   
        
    elif lat_type == "hex":
        ns = range(kpd-5,kpd+6)

    elif lat_type == "trig":
        ns = range(kpd-5,kpd+6)
        
    elif lat_type == "stet":
        ns = range(kpd-5,kpd+6)
        
    elif lat_type == "btet":
        ns = range(kpd-5,kpd+6)
        
    elif lat_type == "so":
        ns = range(kpd-5,kpd+6)
        
    elif lat_type == "co":
        ns = range(kpd-5,kpd+6)
        
    elif lat_type == "bo":
        ns = range(kpd-5,kpd+6)
        
    elif lat_type == "fo":
        ns = range(kpd-5,kpd+6)
        
    elif lat_type == "sm":
        ns = range(kpd-5,kpd+6)
        
    elif lat_type == "so":
        ns = range(kpd-5,kpd+6)
        
    elif lat_type == "tric":
        tmult = 1
        while mult == None:
            test =nt/tmult**3
            if np.allclose(test,100):
                nu = int(ceil(test))
                nl = int(floor(test))
                if (nu-100)<(nl-100):
                    ns = [nu]
                else:
                    ns = [nl]
                mult = tmult
            else:
                tmult += 1
    else:
        raise ValueError("ERROR: unrecognized lattice type: ",lat_type)

    return ns, mult
