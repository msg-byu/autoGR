"""Tests the basis transformaiton over a number of cell sizes for any given cell.
"""

import numpy as np
from opf_python.niggli_lat_id  import niggli_id
from opf_python.pyniggli import reduced_cell
import opf_python.fcc as fcc
import opf_python.sc as sc
import opf_python.bcc as bcc
import opf_python.stet as stet
import opf_python.so as so


def test_cell(U,eps=None):
    """Tests the ability of the niggli basis to preserve the symmetry of
    the parent cell under action of the HNF.

    Args:
        U (numpy array): The parent lattice as columns of vector.
        eps (optional, float): The floating point tolerance, default is 1E-5.

    Returns:
        False if symmetry is not preserved, true if it is. Also returns the 
    """

    B = reduced_cell(U)
    N = B.niggli
    C = B.C

    U_lat_name, U_niggli_n, U_lat_fam = niggli_id(U)

    print("Users niggli number: ", U_niggli_n)
    
    if U_niggli_n == 1:
        for kpd in range(3,23):
            HNFs = fcc.fcc_srHNFs(kpd**3)
            test = _test_case(U,N,C,U_lat_fam,HNFs)
            if not test:
                break

    # if U_niggli_n == 2:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    if U_niggli_n == 3:
        for kpd in range(3,23):
            HNFs = sc.sc_srHNFs(kpd**3)
            test = _test_case(U,N,C,U_lat_fam,HNFs)
            if not test:
                break

    # if U_niggli_n == 4:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    if U_niggli_n == 5:
        for kpd in range(3,23):
            HNFs = bcc.bcc_srHNFs(kpd**3)
            test = _test_case(U,N,C,U_lat_fam,HNFs)
            if not test:
                break
            
    # if U_niggli_n == 6:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 7:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 8:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 9:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 10:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    if U_niggli_n == 11:
        for kpd in range(3,23):
            HNFs = stet.stet_srHNFs(kpd**3)
            test = _test_case(U,N,C,U_lat_fam,HNFs)
            if not test:
                break
            
    # if U_niggli_n == 12:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 13:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 14:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 15:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 16:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 17:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 18:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 19:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 20:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 21:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 22:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 23:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 24:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 25:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 26:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 27:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 28:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 29:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 30:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 31:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    if U_niggli_n == 32:
        for kpd in range(3,23):
            HNFs = so.so_srHNFs(kpd**3)
            test = _test_case(U,N,C,U_lat_fam,HNFs)
            if not test:
                break
    
    # if U_niggli_n == 33:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 34:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 35:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 36:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 37:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 38:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 39:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 40:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 41:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
            
    # if U_niggli_n == 42:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 43:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break
    
    # if U_niggli_n == 44:
    #     for kpd in range(3,23):
    #         HNFs = 
    #         test = _test_case(U,N,C,U_lat_fam,HNFs)
    #         if not test:
    #             break    
    
    return test


def _test_case(U,N,C,crystal_fam,HNFs):
    """Tests a set of HNFs to see if they preserve the symmetry of the
    desired crystal family.

    Args:
        U (numpy array): The users cell as a column matrix.
        N (numpy array): The niggli cell as a column matrix.
        C (numpy array): The transformation from the niggli cell to the users cell. 
        crystal_fam (int): An integer indicating the crystal family of the users cell.
        HNFs (list of lists): A list of the HNF matrices to be testsed. 

    Returns:
        True if the parent cell symmetry is preserved. False otherwise.
    """
    for H in HNFs:
        S = np.dot(np.dot(N,H),np.linalg.inv(C))
        lat_name, niggli_n, lat_fam = niggli_id(S)
        r = np.linalg.inv(np.transpose(U))
        g = np.linalg.inv(np.transpose(S))
        temp = np.round(np.dot(np.linalg.inv(g),r),3)
        if lat_fam > crystal_fam or not np.allclose(temp%1,0):
            print("Failed on HNF: ", H)
            return False

    return True
