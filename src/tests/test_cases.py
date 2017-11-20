"""Tests the srHNF code for each Bravais Lattice."""

import pytest
import numpy as np

def test_sc():
    """Tests the simple cubic lattice srHNF generation.
    """
    from opf_python.sc import sc_3

    srHNFs = []
    for n in range(1,501):
        temp = sc_3(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/sc_1_500_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/sc_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
    
def test_bcc():
    """Tests the body centered cubic lattice srHNF generation.
    """
    from opf_python.bcc import bcc_5

    srHNFs = []
    for n in range(1,501):
        temp = bcc_5(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/bcc_1_500_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/bcc_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
    
def test_fcc():
    """Tests the face centered cubic lattice srHNF generation.
    """
    from opf_python.fcc import fcc_1

    srHNFs = []
    for n in range(1,501):
        temp = fcc_1(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/fcc_1_500_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/fcc_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
        
def test_stet():
    """Tests the simple tetragonal lattice srHNF generation.
    """
    from opf_python.stet import stet_11, stet_21

    # stet 11 check (original)
    srHNFs = []
    for n in range(1,501):
        temp = stet_11(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/stet_1_500_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/stet_11_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 

    # stet 21 check
    srHNFs = []
    for n in range(1,501):
        temp = stet_21(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/stet_21_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute

    # larger cell checks
    for i in range(10):
        k = np.random.random_integers(501,100000)
        print("kpd",k)
        s1 = stet_11(k)
        s2 = stet_21(k)
        assert len(s1)==len(s2)
    
def test_rhom():
    """Tests the trigonal lattice srHNF generation.
    """
    from opf_python.rhom import rhom_9, rhom_24, rhom_4_2

    # niggli 9 (original) check
    srHNFs = []
    for n in range(1,501):
        temp = rhom_9(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/rhom_1_500_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/rhom_9_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute
        
    # niggli 24 check
    srHNFs = []
    for n in range(1,501):
        temp = rhom_24(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/rhom_24_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 

    # niggli 2 and 4 check
    srHNFs = []
    for n in range(1,501):
        temp = rhom_4_2(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/rhom_2_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    for t in srHNFs:
        assert t in brute 
    
    # larger cell checks
    for i in range(10):
        k = np.random.random_integers(501,100000)
        print("kpd",k)
        s1 = rhom_9(k)
        s2 = rhom_4_2(k)
        s3 = rhom_24(k)
        assert len(s1)==len(s2)
        assert len(s1)==len(s3)
    
def test_hex():
    """Tests the hexagonal lattice srHNF generation.
    """
    from opf_python.hx import hex_12, hex_22

    srHNFs = []
    for n in range(1,501):
        temp = hex_12(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/hex_1_500_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/hex_12_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute

    # test niggli cell 22
    srHNFs = []
    for n in range(1,501):
        temp = hex_22(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/hex_22_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute    
    
    # larger cell checks
    for i in range(10):
        k = np.random.random_integers(501,100000)
        print("kpd",k)
        s1 = hex_12(k)
        s2 = hex_22(k)
        assert len(s1)==len(s2)
    
def test_body_tet():
    """Tests the body centered tetragonal lattince srHNF generation.
    """
    from opf_python.body_tet import body_tet_15, body_tet_6, body_tet_7, body_tet_18

    # test body_tet 15 (original)
    srHNFs = []
    for n in range(1,501):
        temp = body_tet_15(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/body_tet_1_500_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/body_tet_15_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
                
    for t in srHNFs:
        assert t in brute

    # test body_tet 6 (original)
    srHNFs = []
    for n in range(1,501):
        temp = body_tet_6(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)
                
    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/body_tet_6_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
                
    for t in srHNFs:
        assert t in brute

    # test body_tet 7 (original)
    srHNFs = []
    for n in range(1,501):
        temp = body_tet_7(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)
                
    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/body_tet_7_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
                
    for t in srHNFs:
        assert t in brute

    # test body_tet 18 (original)
    srHNFs = []
    for n in range(1,501):
        temp = body_tet_18(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)
                
    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/body_tet_18_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
                
    for t in srHNFs:
        assert t in brute

    # larger cell checks
    for i in range(10):
        k = np.random.random_integers(501,100000)
        print("kpd",k)
        s1 = body_tet_15(k)
        s2 = body_tet_18(k)
        s3 = body_tet_6(k)
        s4 = body_tet_7(k)
        assert len(s1)==len(s2)
        assert len(s1)==len(s3)
        assert len(s1)==len(s4)            
        
def test_so():
    """Tests the simple orthorhombic lattice srHNF generation.
    """
    from opf_python.so import so_32

    # test niggli 32 (original)
    srHNFs = []
    for n in range(1,501):
        temp = so_32(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/so_1_500_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/so_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 
    
def test_base_ortho():
    """Tests the base centered orthorhombic lattice srHNF generation.
    """
    from opf_python.base_ortho import base_ortho_38_13, base_ortho_23, base_ortho_40, base_ortho_36

    # niggli cell 38 and 13 (original)
    srHNFs = []
    for n in range(1,501):
        temp = base_ortho_38_13(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/base_ortho_1_500_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_ortho_38_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 
    
    # niggli cell 23
    srHNFs = []
    for n in range(1,501):
        temp = base_ortho_23(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_ortho_23_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 
    
    # niggli cell 40
    srHNFs = []
    for n in range(1,501):
        temp = base_ortho_40(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_ortho_40_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute
    
    # niggli cell 36
    srHNFs = []
    for n in range(1,501):
        temp = base_ortho_36(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_ortho_36_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute

    # larger cell checks
    for i in range(10):
        k = np.random.random_integers(501,100000)
        print("kpd",k)
        s1 = base_ortho_38_13(k)
        s2 = base_ortho_23(k)
        s3 = base_ortho_40(k)
        s4 = base_ortho_36(k)
        assert len(s1)==len(s2)
        assert len(s1)==len(s3)    
        assert len(s1)==len(s4)    
    
def test_body_ortho():
    """Tests the body centered orthorhombic lattice srHNF generation.
    """
    from opf_python.body_ortho import body_ortho_19, body_ortho_8, body_ortho_42

    # niggli cell 19 (original)
    srHNFs = []
    for n in range(1,501):
        temp = body_ortho_19(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/body_ortho_1_500_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/body_ortho_19_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 
    
    # niggli cell 8
    srHNFs = []
    for n in range(1,501):
        temp = body_ortho_8(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/body_ortho_8_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 
    
    # niggli cell 42
    srHNFs = []
    for n in range(1,501):
        temp = body_ortho_42(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/body_ortho_42_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 
    
    # larger cell checks
    for i in range(10):
        k = np.random.random_integers(501,100000)
        print("kpd",k)
        s1 = body_ortho_19(k)
        s2 = body_ortho_8(k)
        s3 = body_ortho_42(k)
        assert len(s1)==len(s2)
        assert len(s1)==len(s3)
    
def test_face_ortho():
    """Tests the face centered orthorhombic lattice srHNF generation.
    """
    from opf_python.face_ortho import face_ortho_26, face_ortho_16

    # niggli cell 26 (original)
    srHNFs = []
    for n in range(1,501):
        temp = face_ortho_26(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/face_ortho_1_500_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/face_ortho_26_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 
    
    # niggli cell 16
    srHNFs = []
    for n in range(1,501):
        temp = face_ortho_16(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/face_ortho_16_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 
    
    # larger cell checks
    for i in range(10):
        k = np.random.random_integers(501,100000)
        print("kpd",k)
        s1 = face_ortho_26(k)
        s2 = face_ortho_16(k)
        assert len(s1)==len(s2)
    
def test_sm():
    """Tests the simple monoclinic lattice srHNF generation.
    """
    from opf_python.sm import sm_33

    srHNFs = []
    for n in range(1,151):
        temp = sm_33(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/sm_1_150_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/sm_1_150_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
    
def test_base_mono():
    """Tests the base centered monoclinic lattice srHNF generation.
    """
    from opf_python.base_mono import base_mono_14, base_mono_27, base_mono_28, base_mono_41, base_mono_43, base_mono_10_17, base_mono_20_25, base_mono_29_30, base_mono_37_39

    # niggli cell 14 (original)
    srHNFs = []
    for n in range(1,201):
        temp = base_mono_14(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/base_mono_1_200_n.out","r") as f:
        n_500 = int(f.readline().strip())

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_mono_14_1_200_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 

    # niggli cell 27
    srHNFs = []
    for n in range(1,201):
        temp = base_mono_27(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_mono_27_1_200_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 

    # niggli cell 28
    srHNFs = []
    for n in range(1,201):
        temp = base_mono_28(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_mono_28_1_200_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 

    # niggli cell 10 and 17
    srHNFs = []
    for n in range(1,201):
        temp = base_mono_10_17(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_mono_10_1_200_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 

    # niggli cell 20 and 25
    srHNFs = []
    for n in range(1,201):
        temp = base_mono_20_25(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_mono_20_1_200_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 

    # niggli cell 29 and 30
    srHNFs = []
    for n in range(1,201):
        temp = base_mono_29_30(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_mono_29_1_200_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 

    # niggli cell 37 and 39
    srHNFs = []
    for n in range(1,201):
        temp = base_mono_37_39(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_mono_37_1_200_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 

    # niggli cell 41
    srHNFs = []
    for n in range(1,201):
        temp = base_mono_41(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_mono_41_1_200_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute 

    # niggli cell 43
    srHNFs = []
    for n in range(1,201):
        temp = base_mono_43(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    assert len(srHNFs) == n_500

    brute = []
    with open("tests/test_output/base_mono_43_1_200_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])
    
    for t in srHNFs:
        assert t in brute
        
    # larger cell checks
    for i in range(10):
        k = np.random.random_integers(201,10000)
        print("kpd",k)
        s1 = base_mono_14(k)
        s2 = base_mono_10_17(k)
        s3 = base_mono_20_25(k)
        s4 = base_mono_27(k)
        s5 = base_mono_28(k)
        s6 = base_mono_29_30(k)
        s7 = base_mono_37_39(k)
        s8 = base_mono_41(k)
        s9 = base_mono_43(k)
        assert len(s1)==len(s2)
        assert len(s1)==len(s3)
        assert len(s1)==len(s4)
        assert len(s1)==len(s5)
        assert len(s1)==len(s6)
        assert len(s1)==len(s7)
        assert len(s1)==len(s8)
        assert len(s1)==len(s9)
