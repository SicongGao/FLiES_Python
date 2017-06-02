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

    def __init__(self):
        self.save_a = 0

    def save(self, a):

        self.save_a = a

        return ERRCODE.SUCCESS
    # calculate the phase function from LUT
    # cb:cb = 1 (overstory),cb = 2 (branch)
    # cb = 3 (stem top/bottom),cb = 4 (understory), cb = 5 (soil surface)
    # cb = 6 (stem side)

    def fpf(self, thi, tho, phr, leaf_reflectance, leaf_transmittance, cb):

        gmrx = np.zeros(3 * 3 * 3, dtype=float).reshape(3, 3, 3)
        gmtx = np.zeros(3 * 3 * 3, dtype=float).reshape(3, 3, 3)
        tr = [0.0] * 3
        tt = [0.0] * 3
        ttt = [0.0] * 3
        trr = [0.0] * 3

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
                    gmrx[l + 1, m + 1, n + 1] = comm.DLT[cb, 1] * comm.G_MRC[i + l, j + m, k + n] + \
                                    comm.DLT[cb, 2] * comm.G_MRB[i + l, j + m, k + n] + \
                                    comm.DLT[cb, 4] * comm.G_MRF[i + l, j + m, k + n]

                    gmtx[l + 1, m + 1, n + 1] = comm.DLT[cb, 1] * comm.G_MTC[i + l, j + m, k + n] + \
                                    comm.DLT[cb, 2] * comm.G_MTB[i + l, j + m, k + n] + \
                                    comm.DLT[cb, 4] * comm.G_MTF[i + l, j + m, k + n]
        # bi-linear over th1 dimension
        for n in range(1, 3):
            for l in range(1, 3):
                tr[l] = gmrx[1, l, n] * (float(i) - th1) + \
                        gmrx[2, l, n] * (th1 - float(i - 1))
                tt[l] = gmtx[1, l, n] * (float(i) - th1) + \
                        gmtx[2, l, n] * (th1 - float(i - 1))

            # bi-linear over th1, th2 dimension
            trr[n] = tr[1] * (float(j) - th2) + tr[2] * (th2 - float(j - 1))
            ttt[n] = tt[1] * (float(j) - th2) + tt[2] * (th2 - float(j - 1))

        # bi-linear over (th1+th2) - ph plan
        gmr = trr[1] * (float(k) - ph) + trr[2] * (ph - float(k - 1))
        gmt = ttt[1] * (float(k) - ph) + ttt[2] * (ph - float(k - 1))

        gm = leaf_reflectance * gmr + leaf_transmittance * gmt
        gfunc = comm.DLT[cb, 1] * comm.GT_BLC[int(th1 * 10.0)] +\
                comm.DLT[cb, 2] * comm.GT_BLB[int(th1 * 10.0)] + \
                comm.DLT[cb, 4] * comm.GT_BLF[int(th1 * 10.0)] + \
                float(comm.DLT[cb, 3] + comm.DLT[cb, 5] + comm.DLT[cb, 6])

        gm = (1.0 / gfunc) * gm / ((leaf_reflectance + leaf_transmittance) * pi)

        return gm

    # cb:cb = 1 (overstory),cb = 2 (branch)
    # cb = 3 (stem top/bottom),cb = 4 (understory), cb = 5 (soil surface)
    # cb = 6 (stem side)
    # a is only used for lambertian reflection from the stem side
    # other case "a" should be 1.0

    # a changed

    def simulate(self, phoCoord, vectCoord, w, leaf_reflectance, leaf_transmittance, cb, a, fd, ichi, ikd, para):
        global count
        if (phoCoord.z <= 1e-5):
            logging.debug("*** count = " + str(count + 1))
        count = 0
        MAX_VALUE = 0.999999
        ua = (sum(comm.U) - comm.U[0]) / comm.N_TS

        ff = [0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0]

        # if x and y is outside area
        # if ((phoCoord.x > comm.X_MAX) or (phoCoord.x < 0) or (phoCoord.y < 0) or (phoCoord.y > comm.Y_MAX)):
        #     return ERRCODE.OUT_OF_RANGE

        objCoord = Position()
        objCoord.x = phoCoord.x - (int(phoCoord.x / comm.X_MAX) - 0.5 + copysign(0.5, phoCoord.x)) * comm.X_MAX
        objCoord.y = phoCoord.y - (int(phoCoord.y / comm.Y_MAX) - 0.5 + copysign(0.5, phoCoord.y)) * comm.Y_MAX
        objCoord.z = phoCoord.z

        # leaf radius (0.1m)
        leafR = 0.1

        vegTrace = VegTrace()
        mc1D = MonteCarlo_1D()

        logging.debug("Vegetation Radiation start...")
        string = "x = " + str(phoCoord.x) + ", y =" + str(phoCoord.y) + ", z =" + str(phoCoord.z)
        logging.debug(string)

        for i in range(1, comm.N_ANG_C + 1): # nangc = nph * nth
            #  preparation of Haple-type hotspot function
            cosa = vectCoord.x * comm.URC_coord[i].x + \
                   vectCoord.y * comm.URC_coord[i].y + \
                   vectCoord.z * comm.URC_coord[i].z
            cosa = copysign(min(cosa, MAX_VALUE), cosa)

            af = acos(cosa)

            th = acos(vectCoord.z)
            ith = int(degrees(th))
            thr = acos(comm.URC_coord[i].z)
            ithr = int(degrees(thr))

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
            # for stem side, the hotspot effect will be ignored so hk=1.0
            ga += comm.DLT[cb, 6] * 1.0e-6
            ch = comm.G_LAI * (comm.DLT[cb, 4] + comm.DLT[cb, 5])
            ch += ua * (comm.DLT[cb, 1] + comm.DLT[cb, 2] + comm.DLT[cb, 3])
            ch += comm.DLT[cb, 6] * 1.0e-6
            ch *= ga * leafR * 0.5

            hk = 1.0 / (1.0 + ch * tan(af * 0.5))
            hk = 1.0 - hk

            vegTrace.trace(objCoord, comm.URC_coord[i])
            sflag = vegTrace.sFlag
            zt = objCoord.z
            objCoord.z = 0.0

            tempCoord = Position()

            taua = mc1D.escape(objCoord, comm.URC_coord[i], 1, ichi, ikd, tempCoord)

            objCoord.z = zt
            thi = acos(vectCoord.z)
            tho = acos(comm.URC_coord[i].z)
            nf = sin(thi) * sin(tho)
            nf = max(nf, 1.0e-8)
            phr = vectCoord.x * comm.URC_coord[i].x + vectCoord.y * comm.URC_coord[i].y
            phr /= nf
            phr = max(-1, phr)
            phr = min(1, phr)
            # print("phr = ", phr)
            phr = acos(phr)

            pf = self.fpf(thi, tho, phr, leaf_reflectance, leaf_transmittance, cb)

            tauc = vegTrace.tau + (- objCoord.z / comm.URC_coord[i].z) * \
                    comm.GT_BLF[ithr] * comm.G_LAI * (comm.DLT[cb, 4] + comm.DLT[cb, 5])
            tauc *= hk

            if (cb == 6):
                rr = sqrt(comm.URC_coord[i].x ** 2 + comm.URC_coord[i].y ** 2)
                rr = max(1.0e-3, rr)
                ph = comm.URC_coord[i].x / rr
                ph = acos(ph)
                a = abs(sin(acos(comm.URC_coord[i].z)) * cos(a - ph))
                self.save(a)

            Id = 0.0
            Id = ff[cb] * w * pf * exp(-tauc) / cos(thr)
            Id += (1.0 - ff[cb]) * w * exp(-tauc) / pi
            Id *= abs(a)
            Id *= (1.0 - fd) * float(sflag)
            # logging.debug("veg rad: id = " + str(Id))
            para.BRF[1, i] += Id
            para.BRF_C[1, i] += Id * comm.DLT[cb, 1]
            para.BRF_S[1, i] += Id * (comm.DLT[cb, 2] + comm.DLT[cb, 3])
            para.BRF_F[1, i] += Id * (comm.DLT[cb, 4] + comm.DLT[cb, 5])

            # Nadir image (nadir lowest point)
            ix = int(objCoord.x * comm.RES) + 1
            iy = int(objCoord.y * comm.RES) + 1

            ix = min(ix, comm.SIZE - 1)
            iY = min(iy, comm.SIZE - 1)

            para.REFL[1, ix, iy] += Id
            para.I_REFL[1, ix, iy] += 1

            Id *= exp(-taua)
            para.BRF[2, i] += Id
            para.BRF_C[2, i] += Id * comm.DLT[cb, 1]
            para.BRF_S[2, i] += Id * comm.DLT[cb, 2] + comm.DLT[cb, 3]
            para.BRF_F[2, i] += Id * comm.DLT[cb, 4] + comm.DLT[cb, 5]

            para.REFL[2, ix, iy] += Id
            para.I_REFL[2, ix, iy] += 1

        self.save(a)
        logging.debug("Vegetation Radiation finish.")
        return ERRCODE.SUCCESS
