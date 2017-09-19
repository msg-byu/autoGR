import numpy as np
import pytest
from opf_python.niggli_lat_id import niggli_id

def test_bo():
    """Tests of a failed body centered orthorhombic case.
    """

    A = [[ 0.07364488,  0.        ,  0.        ],
         [ 0.        ,  0.06377834,  0.        ],
       [ 0.48104638,  0.48104638,  0.96209276]]
    lat_t, nig_n, lat_f = niggli_id(A,1E-10)
    assert lat_t == 'face centered orthorhombic'
    assert lat_f == 5
