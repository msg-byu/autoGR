"""Subtroutines needed by all the systems."""

import numpy as np

def transform_supercells(spHNFs, No, Nu, Co, Cu, O):
    """Finds the symmetry preserving supercells in the users basis.
    
    Args:
        spHNFs (list of lists): The symmetry preserving HNFs.
        No (numpy array): The canonical niggli basis.
        Nu (numpy array): The users niggli basis.
        Co (numpy array): The transformation matrix to the canonical 
            niggli cell.
        Cu (numpy array): The transformation matrix to the users 
            niggli cell.
        O (numpy array): The canonical basis as columns of a matrix.

    Returns:
        spBu (list of lists): The symmetry preserving supercells.
    """
    
    import numpy as np

    spBu = []
    
    for H in spHNFs:
        L = np.dot(np.dot(O,H),np.linalg.inv(O))
        F = np.dot(np.linalg.inv(No),np.dot(np.dot(L,O),Co))
        Bu = np.dot(np.dot(Nu,F),np.linalg.inv(Cu))
        spBu.append(Bu)

    return spBu

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

def find_supercells(U,kpt,exact=False):
    """Generates the symmetry preserving supercells for the users lattice and
    k-point density.

    Args:
        U (numpy.array): The lattice vectors as columns of a matrix.
        kpt (int): The target k-point density desired.
        exact (bool): True if only the target k-point density is to be used.

    Returns:
        supercells (numpy.array): list of symmetry preserving supercells.

    Raises:
        ValueError if the lattice type can't be identified.

    """

    from opf_python.niggli_lat_id import niggli_id
    from opf_python.pyniggli import reduced_cell
    import numpy as np
    
    lat_name, nig_n, lat_fam, basis = niggli_id(U)
    Bu = reduced_cell(U)
    Nu = Bu.niggli
    Cu = Bu.C
    Bo = reduced_cell(basis)
    No = Bu.niggli
    Co = Bu.C
    
    if not exact:
        ns, mult = find_volumes(nig_n,kpt)
    else:
        ns, mult = [kpt], 1.0

    count = 0
    spHNFs = []
    if nig_n == 3:
        from opf_python.sc import sc_3
        for n in ns:
            temp = sc_3(n)
            for t in temp:
                spHNFs.append(t)
        
    elif nig_n == 5:
        from opf_python.bcc import bcc_5
        for n in ns:
            temp = bcc_5(n)
            for t in temp:
                spHNFs.append(t)        
            
    elif nig_n == 1:
        from opf_python.fcc import fcc_1
        for n in ns:
            temp = fcc_1(n)
            for t in temp:
                spHNFs.append(t)        
        
    elif nig_n in [12,22]:
        from opf_python.hx import hex_12, hex_22
        for n in ns:
            if nig_n == 12:
                temp = hex_12(n)
            else:
                temp = hex_22(n)
            if len(temp)> 1 and not exact:
                count +=1
                for t in temp:
                    spHNFs.append(t)
            elif exact:
                for t in temp:
                    spHNFs.append(t)
            if count == 5:
                break
                    
    elif nig_n in [9,4,2,24]:
        from opf_python.rhom import rhom_9, rhom_4_2, rhom_24
        for n in ns:
            if nig_n == 9:
                temp = rhom_9(n)
            elif nig_n == 24:
                temp = rhom_24(n)
            else:
                temp = rhom_4_2(n)                                
            if len(temp)> 1 and not exact:
                count += 1
                for t in temp:
                    spHNFs.append(t)
            elif exact:
                for t in temp:
                    spHNFs.append(t)
            if count == 5:
                break
        
    elif nig_n in [11,21]:
        from opf_python.stet import stet_11, stet_21
        for n in ns:
            if nig_n == 11:
                temp = stet_11(n)
            else:
                temp = stet_21(n)                
            if len(temp)> 1 and not exact:
                count += 1
                for t in temp:
                    spHNFs.append(t)
            elif exact:
                for t in temp:
                    spHNFs.append(t)
            if count == 5:
                break
       
    elif nig_n in [15,7,6,18]:
        from opf_python.body_tet import body_tet_15, body_tet_7, body_tet_6, body_tet_18
        for n in ns:
            if nig_n == 15:
                temp = body_tet_15(n)
            elif nig_n == 18:
                temp = body_tet_18(n)
            elif nig_n == 6:
                temp = body_tet_6(n)
            else:
                temp = body_tet_7(n)
            if len(temp)> 1 and not exact:
                count += 1
                for t in temp:
                    spHNFs.append(t)
            elif exact:
                for t in temp:
                    spHNFs.append(t)
            if count == 5:
                break
        
    elif nig_n == 32:
        from opf_python.so import so_32
        for n in ns:
            temp = so_32(n)
            if len(temp)> 1 and not exact:
                count += 1
                for t in temp:
                    spHNFs.append(t)
            elif exact:
                for t in temp:
                    spHNFs.append(t)
            if count == 5:
                break
        
    elif nig_n in [38,13,23,40,36,]:
        from opf_python.base_ortho import base_ortho_38_13, base_ortho_23, base_ortho_36, base_ortho_40
        for n in ns:
            if nig_n == 23:
                temp = base_ortho_23(n)
            elif nig_n == 36:
                temp = base_ortho_36(n)
            elif nig_n == 40:
                temp = base_ortho_40(n)
            else:
                temp = base_ortho_38_13(n)
            if len(temp)> 1 and not exact:
                count += 1
                for t in temp:
                    spHNFs.append(t)
            elif exact:
                for t in temp:
                    spHNFs.append(t)
            if count == 5:
                break
        
    elif nig_n in [19,8,42]:
        from opf_python.body_ortho import body_ortho_19, body_ortho_8, body_ortho_42
        for n in ns:
            if nig_n == 19:
                temp = body_ortho_19(n)
            elif nig_n == 8:
                temp = body_ortho_8(n)
            else:
                temp = body_ortho_42(n)
            if len(temp)> 1 and not exact:
                count += 1
                for t in temp:
                    spHNFs.append(t)
            elif exact:
                for t in temp:
                    spHNFs.append(t)
            if count == 5:
                break
        
    elif nig_n in [26,16]:
        from opf_python.face_ortho import face_ortho_26, face_ortho_16
        for n in ns:
            if nig_n == 26:
                temp = face_ortho_26(n)
            else:
                temp = face_ortho_16(n)
            if len(temp)> 1 and not exact:
                count += 1
                for t in temp:
                    spHNFs.append(t)
            elif exact:
                for t in temp:
                    spHNFs.append(t)
            if count == 5:
                break
        
    elif nig_n in [33,34,35]:
        from opf_python.sm import sm_33, sm_34, sm_35
        for n in ns:
            if nig_n == 33:
                temp = sm_33(n)
            elif nig_n == 34:
                temp = sm_34(n)
            else:
                temp = sm_35(n)
            if len(temp)> 1 and not exact:
                count += 1
                for t in temp:
                    spHNFs.append(t)
            elif exact:
                for t in temp:
                    spHNFs.append(t)
            if count == 5:
                break
        
    elif nig_n in [14,10,17,20,25,27,28,29,30,41,37,39,43]:
        from opf_python.base_mono import base_mono_14, base_mono_27, base_mono_28, base_mono_41, base_mono_43, base_mono_10_17, base_mono_20_25, base_mono_29_30, base_mono_37_39
        for n in ns:
            if nig_n == 14:
                temp = base_mono_14(n)
            elif nig_n == 27:
                temp = base_mono_27(n)
            elif nig_n == 27:
                temp = base_mono_27(n)
            elif nig_n == 28:
                temp = base_mono_28(n)
            elif nig_n == 41:
                temp = base_mono_41(n)
            elif nig_n == 43:
                temp = base_mono_43(n)
            elif nig_n in [10,17]:
                temp = base_mono_10_17(n)
            elif nig_n in [20,25]:
                temp = base_mono_20_25(n)
            elif nig_n in [29,30]:
                temp = base_mono_29_30(n)
            else:
                temp = base_mono_37_39(n)
            if len(temp)> 1 and not exact:
                count += 1
                for t in temp:
                    spHNFs.append(t)
            elif exact:
                for t in temp:
                    spHNFs.append(t)
            if count == 5:
                break
        
    elif nig_n in [31,44]:
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
                            spHNFs.append(HNF*mult)
        
    else:
        raise ValueError("ERROR: unrecognized lattice type")

    supercells = transform_supercells(spHNFs, No, Nu, Co, Cu, basis)

    return supercells
        
def find_volumes(lat_type,kpd):
    """Finds the allowed n's for a given lattice around the desired k-point density.
    
    Args:
        lat_type (str): The type of lattice.
        kpd (int): the target valume factor.

    Returns:
        ns (list of int): The target valume's to consider.
        mult (int): multiplication factor for triclinic.
    """

    from math import floor, ceil

    mult = None
    ns = []
    
    nt = float(kpd)
    if lat_type == 3 or lat_type == 5:
        nb = int(floor(nt**(1/3))-1)
        nmax = kpd+1000
        while nb**3<nmax:
            nc = nb**3
            for m in [1,2,4]:
                if m*nc>=kpd:
                    ns.append(m*nc)
            nb += 1
        ns.sort()
        ns = ns[:3]
            
    elif lat_type == 1:
        nb = int(floor(nt**(1/3))-1)
        nmax = kpd+1000
        while nb**3<nmax:
            nc = nb**3
            for m in [1,4,16]:
                if m*nc>=kpd:
                    ns.append(m*nc)
            nb += 1
        ns.sort()
        ns = ns[:3]
            
    elif lat_type == 31 or lat_type == 44:
        tmult = 1
        while mult == None:
            test =nt/tmult**3
            if test<=100:
                nl = int(floor(test))
                ns = [nl]
                mult = tmult
            else:
                tmult += 1
    else:
        ns = range(kpd,kpd+11)

    ns.sort()
    return ns, mult
