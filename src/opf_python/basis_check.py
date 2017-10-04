"""Tests the basis transformaiton over a number of cell sizes for any given cell.
"""

import numpy as np
from opf_python.niggli_lat_id  import niggli_id
from opf_python.pyniggli import reduced_cell
from opf_python.universal import transform_supercells
import opf_python.fcc as fcc
import opf_python.sc as sc
import opf_python.bcc as bcc
import opf_python.stet as st
import opf_python.so as so
import opf_python.rhom as rhom
import opf_python.hx as hx
import opf_python.body_tet as bct
import opf_python.base_ortho as cco
import opf_python.body_ortho as bco
import opf_python.face_ortho as fco
import opf_python.sm as sm
import opf_python.base_mono as ccm

def cell_test(U,count,max_kpd_=None):
    """Tests the ability of the niggli basis to preserve the symmetry of
    the parent cell under action of the HNF.

    Args:
        U (numpy array): The parent lattice as columns of vector.
        count (dict): A dictionary to count how many times each niggli 
            cell is used.
        max_kpd_ (optional int): The largest k-point to use.

    Returns:
        False if symmetry is not preserved, true if it is. Also returns the 
            count of the niggli basis used.
    """

    if max_kpd_ is None:
        max_kpd = 100000
    else:
        max_kpd = max_kpd_
    
    temp = reduced_cell(U)
    Nu = temp.niggli
    Cu = temp.C

    U_lat_name, U_niggli_n, U_lat_fam, O = niggli_id(U)
    print("U_niggli_n", U_niggli_n, "U_lat_fam", U_lat_fam)

    tempo = reduced_cell(O)
    No = tempo.niggli
    Co = tempo.C

    count[str(U_niggli_n)] += 1
    
    if U_niggli_n == 1:
        kpds = fcc_kps(max_kpd)
        for kpd in kpds:
            HNFs = fcc.fcc_1(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count

    if U_niggli_n == 2:
        for kpd in range(max_kpd):
            HNFs = rhom.rhom_4_2(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 3:
        kpds = sc_and_bcc_kps(max_kpd)
        for kpd in kpds:
            HNFs = sc.sc_3(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count

    if U_niggli_n == 4:
        for kpd in range(max_kpd):
            HNFs = rhom.rhom_4_2(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 5:
        kpds = sc_and_bcc_kps(max_kpd)
        for kpd in kpds:
            HNFs = bcc.bcc_5(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 6:
        for kpd in range(max_kpd):
            HNFs = bct.body_tet_6(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 7:
        for kpd in range(max_kpd):
            HNFs = bct.body_tet_7(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 8:
        for kpd in range(max_kpd):
            HNFs = bco.body_ortho_8(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 9:
        for kpd in range(max_kpd):
            HNFs = rhom.rhom_9(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 10:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_10_17(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 11:
        for kpd in range(max_kpd):
            HNFs = st.stet_11(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 12:
        for kpd in range(max_kpd):
            HNFs = hx.hex_12(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 13:
        for kpd in range(max_kpd):
            HNFs = cco.base_ortho_38_13(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 14:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_14(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 15:
        for kpd in range(max_kpd):
            HNFs = bct.body_tet_15(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 16:
        for kpd in range(max_kpd):
            HNFs = fco.face_ortho_16(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 17:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_10_17(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 18:
        for kpd in range(max_kpd):
            HNFs = bct.body_tet_18(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 19:
        for kpd in range(max_kpd):
            HNFs = bco.body_ortho_19(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 20:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_20_25(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 21:
        for kpd in range(max_kpd):
            HNFs = st.stet_21(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 22:
        for kpd in range(max_kpd):
            HNFs = hx.hex_22(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 23:
        for kpd in range(max_kpd):
            HNFs = cco.base_ortho_23(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 24:
        for kpd in range(max_kpd):
            HNFs = rhom.rhom_24
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 25:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_20_25(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 26:
        for kpd in range(max_kpd):
            HNFs = fco.face_ortho_26(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 27:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_27(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 28:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_28(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 29:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_29_30(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 30:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_29_30(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 31: # triclinic, nothing to do.
        return True, count
    
    if U_niggli_n == 32:
        for kpd in range(max_kpd):
            HNFs = so.so_32(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 33: # not implemented yet
        return True, count 
        # for kpd in range(max_kpd):
        #     HNFs = sm.sm_33(kpd)
        #     test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
        #     if not test:
        #         return False, count
            
    if U_niggli_n == 34: # not implemented yet
        return True, count 
        # for kpd in range(max_kpd):
        #     HNFs = sm.sm_34(kpd)
        #     test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
        #     if not test:
        #         return False, count
    
    if U_niggli_n == 35: # not implemented yet
        return True, count 
        # for kpd in range(max_kpd):
        #     HNFs = sm.sm_35(kpd)
        #     test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
        #     if not test:
        #         return False, count
            
    if U_niggli_n == 36:
        for kpd in range(max_kpd):
            HNFs = cco.base_ortho_36(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 37:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_37_39(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 38:
        for kpd in range(max_kpd):
            HNFs = cco.base_ortho_38_13(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 39:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_37_39(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 40:
        for kpd in range(max_kpd):
            HNFs = cco.base_ortho_40(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 41:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_41(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
            
    if U_niggli_n == 42:
        for kpd in range(max_kpd):
            HNFs = bco.body_ortho_42(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 43:
        for kpd in range(max_kpd):
            HNFs = ccm.base_mono_43(kpd)
            test = _test_case(U,Nu,Cu,O,No,Co,U_lat_fam,HNFs)
            if not test:
                return False, count
    
    if U_niggli_n == 44: #triclinic
        return True, count
    
    return True, count

def sc_and_bcc_kps(kmax):
    """Lists the k-point densities up to kmax for the sc and bcc cases.

    Args:
        kmax (int): The largest allowed k-point density.

    Returns:
        kpds (list): List of the allowed k-points for sc and bcc systems.
    """

    kpds = []
    for k in range(kmax):
        for m in range(1,2,4):
            p = m*(k**3)
            if p <= kmax and p not in kpds:
                kpds.append(p)

    return kpds

def fcc_kps(kmax):
    """Lists the k-point densities up to kmax for the fcc cases.

    Args:
        kmax (int): The largest allowed k-point density.

    Returns:
        kpds (list): List of the allowed k-points for fcc systems.
    """

    kpds = []
    for k in range(kmax):
        for m in range(1,4,16):
            p = m*(k**3)
            if p <= kmax and p not in kpds:
                kpds.append(p)

    return kpds

def _test_case(U,Nu,Cu,O,No,Co,crystal_fam,HNFs):
    """Tests a set of HNFs to see if they preserve the symmetry of the
    desired crystal family.

    Args:
        U (numpy array): The users cell as a column matrix.
        Nu (numpy array): The users niggli cell as a column matrix.
        Cu (numpy array): The transformation from the users niggli cell to the users cell. 
        No (numpy array): Our niggli cell as a column matrix.
        Co (numpy array): The transformation from our niggli cell to our. 
        O (numpy array): Our cell as a column matrix.
        crystal_fam (int): An integer indicating the crystal family of the users cell.
        HNFs (list of lists): A list of the HNF matrices to be testsed. 

    Returns:
        True if the parent cell symmetry is preserved. False otherwise.
    """
    Bs = transform_supercells(HNFs,No,Nu,Co,Cu,O)
    for i in range(len(Bs)):
        B = Bs[i]
        lat_name, niggli_n, lat_fam, c_b = niggli_id(B)
        r = np.linalg.inv(np.transpose(U))
        g = np.linalg.inv(np.transpose(B))
        temp = np.round(np.dot(np.linalg.inv(g),r),3)
        if lat_fam > crystal_fam or not np.allclose(temp%1,0):
            print('lf',lat_fam,'cf',crystal_fam,'com',temp%1,"HNF",HNFs[i])
            return False

    return True
