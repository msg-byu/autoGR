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
    """Tests the simple cubic lattice srHNF generation.
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
    """Tests the simple cubic lattice srHNF generation.
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
    """Tests the simple cubic lattice srHNF generation.
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
    """Tests the simple cubic lattice srHNF generation.
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
