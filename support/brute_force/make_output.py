import csv
import numpy as np
from niggli_lat_id import niggli_id

As = [np.transpose([[1,0,0],[0,1,0],[0,0,2]]),np.transpose([[0.5,0,0.5],[0,3,0],[0.5,0,-0.5]]),np.transpose([[1.00000000,0.00000000,0.00000000],[-0.50000000,0.86602540,1.63299320],[0.00000000,-1.73205080,1.63299320]]),np.transpose([[1.00000000, 0.00000000, 0.00000000],[-0.50000000, 0.86602540, 0.00000000],[0.00000000, 0.00000000, 3.26598640]]),np.transpose([[1.00000000, 0.00000000, 0.00000000],[-0.50000000, 0.86602540, 1.63299320],[0.00000000, -1.73205080, 1.63299320]]),np.transpose([[1.00000000, 0.00000000, 0.00000000],[0.50000000, -0.86602540, 3.26598640],[0.00000000, -1.73205080, 0.00000000]]),np.transpose([[1.00000000, 0.00000000, 0.00000000],[0.50000000, 4.33012700, 0.00000000],[0.00000000, 0.00000000, 1.63299320]]),np.transpose([[1.00000000, 0.00000000, 0.00000000],[0.50000000, -0.86602540, 4.89897960],[0.00000000, -1.73205080, 0.00000000]]),np.transpose([[0.00000000, -1.73205080, 1.63299320],[0.50000000, 2.59807620, 3.26598640],[1.00000000, 0.00000000, 0.00000000]]),np.transpose([[0.5, 0.5, -0.5],[-0.5, 0.5, 0.5],[1.0, 0.0, 1.0]]),np.transpose([[1.00000000, 0.00000000, 0.00000000],[0.00000000, 0.00000000, 1.00000000],[0.50000000, -1.50000000, 0.50000000]]),np.transpose([[.05, 2.7, 3.3],[0.1, 0.7, 4.5],[.99, .3, 5.4]]),np.transpose([[1, -.1, 0],[-0.3, 1, .3],[-.3, -0.1, -1.5 ]]),np.array([[ 0.03682244,  0.        ,  0.        ],[ 0.        ,  0.06377834,  0.        ],[ 0.96209276,  0.96209276,  1.92418552]]),np.array([[ 0.07364488,  0.        ,  0.        ],[ 0.        ,  0.12755668,  0.        ],[ 0.24052319,  0.24052319,  0.48104638]]),np.array([[ 1.20559446,  0.        ,  0.        ],[ 0.        ,  2.08815085,  0.        ],[ 3.93745511,  3.93745511,  7.87491022]])]

eps = 1E-5

for i  in range(len(As)):

    U = As[i]
    lat_name, nig_n, lat_fam, O = niggli_id(U)
    Bu = reduced_cell(U)
    Nu = Bu.niggli
    Cu = Bu.C
    Bo = reduced_cell(O)
    No = Bo.niggli
    Co = Bo.C

    in_f = open("U.in.{}".format(i+1),"w+")
    w1 = csv.writer(in_f,delimiter="\t")
    for d in U:
        w1.writerow(d)
    in_f.close()

    o_f = open("O.out.{}".format(i+1),"w+")
    w2 = csv.writer(in_f,delimiter="\t")
    for d in O:
        w2.writerow(d)
    o_f.close()
    o_f = open("No.out.{}".format(i+1),"w+")
    w2 = csv.writer(in_f,delimiter="\t")
    for d in No:
        w2.writerow(d)
    o_f.close()
    o_f = open("Co.out.{}".format(i+1),"w+")
    w2 = csv.writer(in_f,delimiter="\t")
    for d in Co:
        w2.writerow(d)
    o_f.close()
    o_f = open("Nu.out.{}".format(i+1),"w+")
    w2 = csv.writer(in_f,delimiter="\t")
    for d in Nu:
        w2.writerow(d)
    o_f.close()
    o_f = open("No.out.{}".format(i+1),"w+")
    w2 = csv.writer(in_f,delimiter="\t")
    for d in No:
        w2.writerow(d)
    o_f.close()

    o_f = open("id.out.{}".format(i+1),"w+")
    o_f.write(str(nig_n))
    o_f.close()
