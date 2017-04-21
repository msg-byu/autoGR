"""Tests the srHNF code for each Bravais Lattice."""

import pytest
import sys

def test_sc():
    """Tests the simple cubic lattice srHNF generation.
    """
    from opf_python.sc import sc_srHNFs

    srHNFs = []
    for n in range(1,501):
        temp = sc_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/sc_1_500_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

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
    """Tests the base centered cubic lattice srHNF generation.
    """
    from opf_python.bcc import bcc_srHNFs

    srHNFs = []
    for n in range(1,501):
        temp = bcc_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/bcc_1_500_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

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
    from opf_python.fcc import fcc_srHNFs

    srHNFs = []
    for n in range(1,501):
        temp = fcc_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/fcc_1_500_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

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
    from opf_python.stet import stet_srHNFs

    srHNFs = []
    for n in range(1,501):
        temp = stet_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/stet_1_500_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

    brute = []
    with open("tests/test_output/stet_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
    
def test_trig():
    """Tests the trigonal lattice srHNF generation.
    """
    from opf_python.trig import trig_srHNFs

    srHNFs = []
    for n in range(1,501):
        temp = trig_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/trig_1_500_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

    brute = []
    with open("tests/test_output/trig_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
    
def test_hex():
    """Tests the hexagonal lattice srHNF generation.
    """
    from opf_python.hx import hex_srHNFs

    srHNFs = []
    for n in range(1,501):
        temp = hex_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/hex_1_500_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

    brute = []
    with open("tests/test_output/hex_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
    
def test_body_tet():
    """Tests the body centered tetragonal lattince srHNF generation.
    """
    from opf_python.body_tet import body_tet_srHNFs

    srHNFs = []
    for n in range(1,501):
        temp = body_tet_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/body_tet_1_500_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

    brute = []
    with open("tests/test_output/body_tet_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
    
def test_so():
    """Tests the simple orthorhombic lattice srHNF generation.
    """
    from opf_python.so import so_srHNFs

    srHNFs = []
    for n in range(1,501):
        temp = so_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/so_1_500_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

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
    from opf_python.base_ortho import base_ortho_srHNFs

    srHNFs = []
    for n in range(1,501):
        temp = base_ortho_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/base_ortho_1_500_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

    brute = []
    with open("tests/test_output/base_ortho_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
    
def test_body_ortho():
    """Tests the body centered orthorhombic lattice srHNF generation.
    """
    from opf_python.body_ortho import body_ortho_srHNFs

    srHNFs = []
    for n in range(1,501):
        temp = body_ortho_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/body_ortho_1_500_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

    brute = []
    with open("tests/test_output/body_ortho_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
    
def test_face_ortho():
    """Tests the face centered orthorhombic lattice srHNF generation.
    """
    from opf_python.face_ortho import face_ortho_srHNFs

    srHNFs = []
    for n in range(1,501):
        temp = face_ortho_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/face_ortho_1_500_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

    brute = []
    with open("tests/test_output/face_ortho_1_500_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
    
def test_sm():
    """Tests the simple monoclinic lattice srHNF generation.
    """
    from opf_python.sm import sm_srHNFs

    srHNFs = []
    for n in range(1,151):
        temp = sm_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/sm_1_150_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

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
    from opf_python.base_mono import base_mono_srHNFs

    srHNFs = []
    for n in range(1,201):
        temp = base_mono_srHNFs(n)
        for t in temp:
            if len(t) >0:
                srHNFs.append(t)

    with open("tests/test_output/base_mono_1_200_n.out","r") as f:
        n = int(f.readline().strip())

    assert len(srHNFs) == n

    brute = []
    with open("tests/test_output/base_mono_1_200_srHNFs.out","r") as f:
        HNF = []
        for line in f:
            if len(line.strip().split()) == 0:
                brute.append(HNF)
                HNF = []
            else:
                HNF.append([int(i) for i in line.strip().split()])

    
    for t in srHNFs:
        assert t in brute 
