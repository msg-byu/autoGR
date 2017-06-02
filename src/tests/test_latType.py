"""Tests the srHNF code for each Bravais Lattice."""

import pytest
import sys
from opf_python.lat_type import lat_type
import numpy as np

def test_sc():
    """Tests identification of simple cubic lattices.
    """

    As = [[[1,0,0],[0,1,0],[0,0,1]], [[1,0,0],[1,1,0],[0,0,1]],
          [[1,0,1],[1,1,0],[0,0,1]],[[1,1,1],[0,1,1],[1,0,1]],[[1,0,2],[0,1,0],[0,1,1]]]

    for A in As:
        name, basis = lat_type(A)

        assert name=="sc"

def test_bcc():
    """Tests identification of body centered cubic lattices.
    """

    As = [[[-1,1,1],[1,-1,1],[1,1,-1]],[[1,0,0],[0,1,0],[0.5,0.5,0.5]],
          [[1,0,0],[0.5,0.5,0.5],[0,0,1]],[[0.5,-0.5,0.5],[0,1,0],[0,0,1]],
          [[0,0,2],[1,-1,1],[1,1,-1]],[[-1,1,1],[2,0,0],[0,2,0]],
          [[1,-1,3],[1,-1,1],[1,1,-1]]]

    for A in As:
        name, basis = lat_type(A)

        assert name=="bcc"

def test_fcc():
    """Tests the identification of face centered cubic lattices.
    """

    As = [[[0,1,1],[1,0,1],[1,1,0]],[[1,0,0],[0.5,0,0.5],[0.5,0.5,0]],
          [[0,0,1],[0.5,0,0.5],[0.5,0.5,0]],[[0,0.5,0.5],[0,1,0],[0.5,0.5,0]],
          [[0,0.5,0.5],[0.5,0,0.5],[1,0,0]],[[0,0.5,0.5],[0.5,0,0.5],[0,0,1]],
          [[1,1,2],[1,0,1],[2,1,1]],[[2,2,2],[1,1,2],[1,2,1]],
          [[0,3.82669513e+00,3.82669513e+00],[3.82669513e+00,0,3.82669513e+00],
           [3.82669513e+00,3.82669513e+00,0]]]

    for A in As:
        name, basis = lat_type(A)

        assert name=="fcc"

def test_hex():
    """Tests the identification of the hexagonal lattices.
    """

    a1 = np.array([1,0,0])
    a2 = np.array([0.5, -0.8660254037844386, 0])
    a3 = np.array([0,0,2])
    As = [[a1,[-0.5, -0.8660254037844386, 0],a3],[a1,a2,a3],[a1+a2,a2,a1+a3],
          [a1+a2+a3,a2+a3,a1+a3],[a1+a2+a2,a2,a3+a2],[a1,[-0.5, -0.866025, 0],a3]]

    for A in As:
        name, basis = lat_type(A)

        assert name == "hex"

def test_trig():
    """Tests the identification of trigonal lattices.
    """

    a1 = np.array([1,2,2])
    a2 = np.array([2,1,2])
    a3 = np.array([2,2,1])

    As = [[a1,a2,a3],[a2+a1,a2+a3,a1+a3],[-a1,a2+a1+a1,a3-a1],[a1+a2+a3,a1+a2,a2+a3],
          [a1+a1,a2+a2,a3+a3]]

    for A in As:
        name, basis = lat_type(A)

        assert name == "trig"

def test_simple_tet():
    """Tests the identification of the simple tetragonal lattices.
    """

    a1 = np.array([1,0,0])
    a2 = np.array([0,1,0])
    a3 = np.array([0,0,2])

    As = [[a1,a2,a3],[a1+a2,a2,a3+a1],[a1+a2+a3,a2+a3,a1+a2],[a1,a2+a1,a3+a3],[-a1,-a2+a3,a3]]

    for A in As:
        name, basis = lat_type(A)

        assert name == "stet"

def test_btet():
    """Tests the identfication of the body centered tetragonal lattices.
    """

    a1 = np.array([-1,1,2])
    a2 = np.array([1,-1,2])
    a3 = np.array([1,1,-2])

    As = [[a1,a2,a3],[[2,0,0],[0,2,0],a1],[[2,0,0],a1,[0,0,4]],[a1,[0,2,0],[0,0,4]],
          [a1+a2,a2,a1+a3],[a1+a2+a3,-a2,a3],[a1,a2+a3,a1+a3]]

    for A in As:
        name, basis = lat_type(A)

        assert name == "btet"

def test_so():
    """Tests the identification of simple orthorhombic lattices.
    """

    a1 = np.array([1,0,0])
    a2 = np.array([0,2,0])
    a3 = np.array([0,0,3])

    As = [[a1,a2,a3],[a1+a2,a2,a3],[-a1,a2+a3,-a3],[a1+a2+a3,a2,a3],[a1,a2+a1+a2,a3]]

    for A in As:
        name, basis = lat_type(A)

        assert name == "so"

def test_co():
    """Tests the identification of base centered orthorhombic lattices.
    """

    a1 = np.array([0.5,1,0])
    a2 = np.array([0.5,-1,0])
    a3 = np.array([0,0,3])

    As = [[a1,a2,a3],[a1,[0,2,0],[0,0,3]],[[1,0,0],a1,[0,0,3]],[a1+a2,a1+a3,a3],[a1+a2+a3,a2,a3],
          [a1+a2+a2,a2+a3,a3],[a1+a3,a2,a2+a3]]

    for A in As:
        name, basis = lat_type(A)

        assert name == "co"

def test_bo():
    """Tests the identification of body centeretd orthorhombic lattices.
    """

    a1 = np.array([-0.5,1,1.5])
    a2 = np.array([0.5,-1,1.5])
    a3 = np.array([0.5,1,-1.5])

    As = [[a1,a2,a3],[[1,0,0],[0,2,0],a1],[[1,0,0],a1,[0,0,3]],[a1,[0,2,0],[0,0,3]],[a1+a2,a2,a3]
          ,[a1+a2+a3,a2,a3],[a1+a3,a2+a3,a3],[a1+a2+a2,a2+a3,a3+a1]]

    for A in As:
        name, basis = lat_type(A)

        assert name == "bo"

def test_fo():
    """Tests the identification of the face centered orthorhombic lattices.
    """

    a1 = np.array([0,1,1.5])
    a2 = np.array([0.5,0,1.5])
    a3 = np.array([0.5,1,0])

    As = [[a1,a2,a3],[[1,0,0],a2,a3],[[0,0,3],a2,a3],[a1,[0,2,0],a3],[a1,a2,[1,0,0]],
          [a1,a2,[0,0,3]],[a1+a2,a2,a3],[a1+a2+a3,a2+a3,a3],[-a1,a2-a1,a3+a1],[a1,a2+a1,a2+a3]]

    for A in As:
        name, basis = lat_type(A)

        assert name == "fo"

def test_sm():
    """Tests the identification of the simple monoclinc lattices.
    """

    a1 = np.array([2,0,0])
    a2 = np.array([0,2,0])
    a3 = np.array([0.5,0,2])

    As = [[a1,a2,a3],[a1+a2+a3,a2,a3],[a1+a2,a2,a2+a3],[a1,a2+a3,a1+a3],[a1+a2+a3,a2-a3,a3],
          [a1,a2,a3+a2-a1]]

    for A in As:
        name, basis = lat_type(A)

        assert name == "sm"

def test_cm():
    """Tests the identification of the base centered monoclinic lattices.
    """

    a1 = np.array([1,1,0])
    a2 = np.array([1,-1,0])
    a3 = np.array([0.5,0,2])

    As = [[a1,a2,a3],[a1,[0,2,0],a3],[[2,0,0],a2,a3],[a1+a2+a3,a2,a3],[a1+a3,a1+a2,a2+a3],
          [a1-a2,a2,a3-a1],[a1,a2+a3,a3],[a1+a2+a2,a1+a3,a2]]

    for A in As:
        name, basis = lat_type(A)

        assert name == "cm"
        

def test_tric():
    """Tests the identification of the triclinic lattices.
    """

    As = [[[1,0,1],[1,2.5,0],[1.5,2,3]],[[1,2,0],[2,7,4],[1,2,3]],[[0.5,0.5,0.5],[1,0,7],[0.5,2,3]],
          [[.1,2,2],[2,.3,2],[2,2,.5]],[[2,9,2],[1,3,5],[.1,4,11]],[[-2,3,4],[2,-3,1],[1,3,-1]]]

    for A in As:
        name, basis = lat_type(A)

        assert name == "tric"
