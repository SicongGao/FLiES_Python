from math import *
import ERRCODE
from G_Function import G_Function
import common as comm
import numpy as np


def ipf():

    th1 = th2 = ph = 0.0
    gmr = gmt = 0.0

    rfile = ["..\data\gmr_uni.txt", "..\data\gmr_plano.txt", "..\data\gmr_erect.txt"]
    tfile = ["..\data\gmt_uni.txt", "..\data\gmt_plano.txt", "..\data\gmt_erect.txt"]

    # leaf angle distribution
    # mc = 1
    # mb = 1
    # mf = 1
    # op = 1  kernel calculation, else read from txt file
    op = 2

    # generate by itself
    if (op == 1):
        for i in range(1, 19):
            print("i = ", i)
            th1 = radians(float(i - 1) * 10.0)

            for j in range(1, 19):
                print("j = ", j)
                th2 = radians(float(j - 1) * 10.0)

                for k in range(1, 37):
                    print("k = ", k)
                    ph = radians(float(k - 1) * 10.0)

                    # canopy
                    gmkernel(gmr, gmt, th1, th2, ph, comm.M_C)
                    comm.G_MRC[i, j, k] = gmr
                    comm.G_MTC[i, j, k] = gmt

                    # branch
                    gmkernel(gmr, gmt, th1, th2, ph, comm.M_B)
                    comm.G_MRB[i, j, k] = gmr
                    comm.G_MTB[i, j, k] = gmt

                    # forest floor
                    gmkernel(gmr, gmt, th1, th2, ph, comm.M_F)
                    comm.G_MRF[i, j, k] = gmr
                    comm.G_MTF[i, j, k] = gmt

        # for reflection
        # write rfiles
        f = open(rfile[0], "w")
        f.write(comm.G_MRC)
        f.close()

        f = open(rfile[1], "w")
        f.write(comm.G_MRB)
        f.close()

        f = open(rfile[2], "w")
        f.write(comm.G_MRF)
        f.close()

        # for transmission
        # write tfiles
        f = open(tfile[0], "w")
        f.write(comm.G_MTC)
        f.close()

        f = open(tfile[1], "w")
        f.write(comm.G_MTB)
        f.close()

        f = open(tfile[2], "w")
        f.write(comm.G_MTF)
        f.close()
    # read from files
    else:
        # for reflection
        contentFile = np.loadtxt(rfile[0])
        comm.G_MRC = transferToMatrix(contentFile)

        contentFile = np.loadtxt(rfile[1])
        comm.G_MRB = transferToMatrix(contentFile)

        contentFile = np.loadtxt(rfile[2])
        comm.G_MRF = transferToMatrix(contentFile)

        # for transmission
        contentFile = np.loadtxt(tfile[0])
        comm.G_MTC = transferToMatrix(contentFile)

        contentFile = np.loadtxt(tfile[1])
        comm.G_MTB = transferToMatrix(contentFile)

        contentFile = np.loadtxt(tfile[2])
        comm.G_MTF = transferToMatrix(contentFile)

    return ERRCODE.SUCCESS


# Mathmatics is summarized in
# Shultis and Myneni (1988)
# J. Q. Spec. Radi. Trans., 39(2), p. 115-129
def gmkernel(gmr, gmt, th1, th2, ph, m):

    gmr = 0.0
    gmt = 0.0
    phl = 0.0   # phl is radian

    gFunction = G_Function()

    gmrp = [0.0] * 362
    gmtp = [0.0] * 362
    gmrt = [0.0] * 182
    gmtt = [0.0] * 182

    for i in range(1, 361):
        phl = radians(float(i) - 1.0)

        for j in range(1, 91):
            thl1 = radians(float(j) - 1.0)
            af1 = faf(th1, thl1, phl)
            af2 = faf(th2, thl1, abs(phl - ph))
            cosf = cos(af1) * cos(af2)

            if (cosf < 0.0):
                gmrt[j] = abs(cosf)
                gmtt[j] = 0.0
            else:
                gmrt[j] = 0.0
                gmtt[j] = cosf

        # integration over theta
        for j in range(1, 90):
            thl1 = radians(float(j) - 1)
            gl1 = gFunction.fgl(thl1, m) / sin(max(thl1, 1.0e-8))
            thl2 = thl1 + pi / 180
            gl2 = gFunction.fgl(thl2, m) / sin(max(thl2, 1.0e-8))

            gmrp[i] += (gl2 * gmrt[j + 1] * sin(thl2) + gl1 * gmrt[j] * sin(thl1)) * pi / (2.0 * 2.0 * 90.0)
            gmtp[i] += (gl2 * gmtt[j + 1] * sin(thl2) + gl1 * gmtt[j] * sin(thl1)) * pi / (2.0 * 2.0 * 90.0)

    for i in range(1, 360):
        gmr += (gmrp[i + 1] + gmrp[i]) * 2.0 * pi / (2.0 * 360.0)
        gmr += (gmtp[i + 1] + gmtp[i]) * 2.0 * pi / (2.0 * 360.0)

    gmr /= (2.0 * pi)
    gmt /= (2.0 * pi)

    return ERRCODE.SUCCESS


# calculation of alpha
def faf(th1, th2, ph):

    data = sin(th1) * cos(th2) + sin(th1) * sin(th2) * cos(ph)
    data = copysign(min(abs(data), 1.0 - 1.0e-10), data)
    data = acos(data)

    return data


def transferToMatrix(origin):

    destin = np.zeros(20 * 20 * 38, dtype=float).reshape(20, 20, 38)
    for i in range(0, 19):
        for j in range(0, 19):
            for k in range(0, 37):
                destin[i + 1][j + 1][k + 1] = origin[i*19+j][k]

    return destin
