"""Tests the srHNF code for each Bravais Lattice."""

# import pytest
# import sys
# from opf_python.universal import find_srBs
# import numpy as np

# def test_sc():
#     """Tests the simple cubic lattice srHNF generation.
#     """

#     A = np.transpose([[1,0,0],[0,1,1],[0,0,1]])
#     Bs, spHNFs = find_srBs(A,1000)

#     brute = []
#     with open("tests/test_output/sc_A1_1_2000_Bs.out","r") as f:
#         B = []
#         for line in f:
#             if len(line.strip().split()) ==0:
#                 brute.append(B)
#                 B = []
#             else:
#                 B.append([float(i) for i in line.strip().split()])

#     In = []
#     for t in Bs:
#         for j in brute:
#             temp = np.dot(t,np.linalg.inv(j))
#             if np.allclose(abs(np.linalg.det(temp)),1) and np.allclose(temp%1,0):
#                 In.append(True)
#                 break
        

#     assert len(In) == len(Bs)
            
