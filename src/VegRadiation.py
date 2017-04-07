import ERRCODE
import common as comm
from Position import Position
from math import *
from VegTrace import VegTrace
from MonteCarlo_1D import MonteCarlo_1D
import numpy as np
import logging

count = 0
class VegRadiation:

    save_a = 0


    def __init__(self):
        self.save_a = 0

    def save(self, a):
        global count
        self.save_a = a
        count += 1
        logging.debug("count = " + str(count))
        return ERRCODE.SUCCESS
    # calculate the phase function from LUT
    # cb:cb = 1 (overstory),cb = 2 (branch)
    # cb = 3 (stem top/bottom),cb = 4 (understory), cb = 5 (soil surface)
    # cb = 6 (stem side)

    def fpf(self, thi, tho, phr, lr, lt, cb):

        gmrx = np.zeros(2 * 2 * 2, dtype=float).reshape(2, 2, 2)
        gmtx = np.zeros(2 * 2 * 2, dtype=float).reshape(2, 2, 2)
        tr = [0.0] * 2
        tt = [0.0] * 2
        ttt = [0.0] * 2
        trr = [0.0] * 2

        # TODO check the math.degree or math.radians function.
        # change angles rad to degree
        th1 = min(degrees(thi), 179.999999)
        th2 = min(degrees(tho), 179.999999)
        ph = degrees(phr)

        # all angles are divided by 10.
        th1 /= 10
        th2 /= 10
        ph /= 10

        i = int(trunc(th1)) + 1
        j = int(trunc(th2)) + 1
        k = int(trunc(ph)) + 1

        for l in range(2):
            for m in range(2):
                for n in range(2):
                    gmrx[l, m, n] = comm.DLT[cb, 1] * comm.G_MRC[i + l + 1, j + m + 1, k + n + 1] + \
                                    comm.DLT[cb, 2] * comm.G_MRB[i + l + 1, j + m + 1, k + n + 1] + \
                                    comm.DLT[cb, 4] * comm.G_MRF[i + l + 1, j + m + 1, k + n + 1]

                    gmtx[l, m, n] = comm.DLT[cb, 1] * comm.G_MRC[i + l + 1, j + m + 1, k + n + 1] + \
                                    comm.DLT[cb, 2] * comm.G_MRB[i + l + 1, j + m + 1, k + n + 1] + \
                                    comm.DLT[cb, 4] * comm.G_MRF[i + l + 1, j + m + 1, k + n + 1]

        # bi-linear over th1 dimension
        for n in range(2):
            for l in range(2):
                tr[l] = gmrx[0, l, n] * (float(i) - th1) + \
                        gmrx[1, l, n] * (th1 - float(i - 1))
                tt[l] = gmrx[0, l, n] * (float(i) - th1) + \
                        gmrx[1, l, n] * (th1 - float(i - 1))

            # bi-linear over th1, th2 dimension
            trr[n] = tr[0] * (float(j) - th2) + tr[1] * (th2 - float(j-1))
            ttt[n] = tt[0] * (float(j) - th2) + tt[1] * (th2 - float(j - 1))

        # bi-linear over (th1+th2) - ph plan
        gmr = trr[0] * (float(k) - ph) + trr[1] * (ph - float(k - 1))
        gmt = ttt[0] * (float(k) - ph) + trr[1] * (ph - float(k - 1))

        gm = lr * gmr + lt * gmt
        gfunc = comm.DLT[cb, 1] * comm.GT_BLC[int(th1 * 10.0)] +\
                comm.DLT[cb, 2] * comm.GT_BLB[int(th1 * 10.0)] + \
                comm.DLT[cb, 4] * comm.GT_BLF[int(th1 * 10.0)] + \
                float(comm.DLT[cb, 3] + comm.DLT[cb, 5] + comm.DLT[cb, 6])

        gm = (1.0 / gfunc) * gm / ((lr + lt) * pi)

        return gm

    # cb:cb = 1 (overstory),cb = 2 (branch)
    # cb = 3 (stem top/bottom),cb = 4 (understory), cb = 5 (soil surface)
    # cb = 6 (stem side)
    # a is only used for lambertian reflection from the stem side
    # other case "a" should be 1.0

    # a changed

    def simulate(self, phoCoord, vectCoord, w, lr, lt, cb, a, fd, ichi, ikd):

        if (phoCoord.z <= 1e-5):
            logging.debug("***count = " + str(count))

        MAX_VALUE = 0.999999
        ua = sum(comm.U) / comm.N_TS

        ff = (1.0, 1.0, 0.0, 1.0, 0.0, 0.0)

        # if x and y is outside area
        objCoord = Position()
        objCoord.x = phoCoord.x - (int(phoCoord.x / comm.X_MAX) - 0.5 + copysign(0.5, phoCoord.x)) * comm.X_MAX
        objCoord.y = phoCoord.y - (int(phoCoord.y / comm.Y_MAX) - 0.5 + copysign(0.5, phoCoord.y)) * comm.Y_MAX
        objCoord.z = phoCoord.z

        # leaf radius (0.1m)
        leafR = 0.1

        vegTrace = VegTrace()
        mc1D = MonteCarlo_1D()

        logging.debug("Vegetation Radiation start...")

        for i in range(comm.N_ANG_C): # nangc = nph * nth
            #  preparation of Haple-type hotspot function
            cosa = vectCoord.x * comm.URC_coord[i].x + \
                   vectCoord.y * comm.URC_coord[i].y + \
                   vectCoord.z * comm.URC_coord[i].z
            cosa = copysign(min(cosa, MAX_VALUE), cosa)

            af = acos(cosa)

            th = acos(vectCoord.z)
            ith = int(radians(th))
            thr = acos(comm.URC_coord[i].z)
            ithr = int(radians(thr))

            # Hapke, types hot spot function
            # af is converted to the opposite angle of scattering angle
            # See Kobayashi and Iwabuchi (2008) Remote Sensing of Environment
            af = pi - af
            ga = comm.DLT[cb, 1] * (comm.GT_BLC[ith] + comm.GT_BLC[ithr]) * 0.5
            ga += comm.DLT[cb, 2] * (comm.GT_BLB[ith] + comm.GT_BLB[ithr]) * 0.5
            ga += comm.DLT[cb, 3] * (comm.GT_BLC[ith] + comm.GT_BLC[ithr]) * 0.25
            ga += comm.DLT[cb, 3] * (comm.GT_BLB[ith] + comm.GT_BLB[ithr]) * 0.25
            ga += comm.DLT[cb, 4] * (comm.GT_BLF[ith] + comm.GT_BLF[ithr]) * 0.5
            ga += comm.DLT[cb, 5] * (comm.GT_BLF[ith] + comm.GT_BLF[ithr]) * 0.5

            ch = comm.G_LAI * (comm.DLT[cb, 4] + comm.DLT[cb, 5])
            ch += ua * (comm.DLT[cb, 1] + comm.DLT[cb, 2] + comm.DLT[cb, 3])
            ch *= ga * leafR * 0.5
            ch = 1.0 / ch
            hk = 1.0 / (1.0 + ch * tan(af * 0.5))
            hk = 1.0 - hk

            vegTrace.trace(objCoord, comm.URC_coord[i])
            sflag = vegTrace.sFlag
            zt = objCoord.z
            objCoord.z = 0.0

            tempCoord = Position()
            # TODO tempCoord doesn't used
            taua = mc1D.escape(objCoord, comm.URC_coord[i], 1, ichi, ikd, tempCoord)

            thi = acos(vectCoord.z)
            tho = acos(comm.URC_coord[i].z)
            nf = sin(thi) * sin(tho)
            nf = max(nf, 1.0e-8)
            phr = vectCoord.x * comm.URC_coord[i].x + vectCoord.y * comm.URC_coord[i].y
            phr /= nf
            phr = acos(phr)

            pf = self.fpf(thi, tho, phr, lr, lt, cb)

            tauc = vegTrace.tau + (- objCoord.z / comm.URC_coord[i].z) * \
                    comm.GT_BLF[ithr] * comm.G_LAI * (comm.DLT[cb, 4] + comm.DLT[cb, 5])
            tauc *= hk

            if (cb == 6):
                rr = sqrt(comm.URC_coord[i].x ** 2 + comm.URC_coord[i].y ** 2)
                rr = max(1.0e-3, rr)
                ph = comm.URC_coord[i].x / rr
                ph = acos(ph)
                a = abs(sin(acos(comm.URC_coord[i].z)) * cos(a - ph))

            Id = 0.0
            Id = ff[cb] * w * pf * exp(-tauc) / cos(thr)
            Id += (1.0 - ff[cb]) * w * exp(-tauc) / pi
            Id *= abs(a)
            Id *= (1.0 - fd) * float(sflag)

            comm.BRF[1, i] += Id
            comm.BRF_C[1, i] += Id * comm.DLT[cb, 1]
            comm.BRF_S[1, i] += Id * comm.DLT[cb, 2] + comm.DLT[cb, 3]
            comm.BRF_F[1, i] += Id * comm.DLT[cb, 4] + comm.DLT[cb, 5]

            # Nadir image (nadir lowest point)
            ix = int(objCoord.x * comm.RES) + 1
            iy = int(objCoord.y * comm.RES) + 1

            ix = min(ix, comm.SIZE - 1)
            iY = min(iy, comm.SIZE - 1)

            comm.REFL[1, ix, iy] += Id
            comm.I_REFL[1, ix, iy] += 1

            Id *= exp(-taua)
            comm.BRF[2, i] += Id
            comm.BRF_C[2, i] += Id * comm.DLT[cb, 1]
            comm.BRF_S[2, i] += Id * comm.DLT[cb, 2] + comm.DLT[cb, 3]
            comm.BRF_F[2, i] += Id * comm.DLT[cb, 4] + comm.DLT[cb, 5]

            comm.REFL[2, ix, iy] += Id
            comm.I_REFL[2, ix, iy] += 1

        self.save(a)
        global count
        logging.debug("Vegetation Radiation finish.")
        return ERRCODE.SUCCESS
