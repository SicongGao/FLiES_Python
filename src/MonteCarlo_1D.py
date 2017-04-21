import ERRCODE
import common as comm
from Position import Position
from math import *
from VegTrace import VegTrace
from Random import Random
import logging


class MonteCarlo_1D:

    weight = -1
    tau = -1

    def __init__(self):
        self.weight = -1
        self.tau = -1
        self.g0fwd = 0.0
        self.g0bwd = 0.0
        self.g1 = 0.0
        self.g1fwd = 0.0
        self.g1bwd = 0.0
        self.g2 = 0.0
        self.g2fwd = 0.0
        self.g2bwd = 0.0
        self.fd = 0.0
        self.ftf = 0.0
        self.ftb = 0.0
        self.g1t = 0.0
        self.g2t = 0.0
        self.angtf = 0.0
        self.angtb = 0.0
        self.iangtf = 0.0
        self.iangtb = 0.0
        self.random_r2 = 0.0
        self.random_r = 0.0
        self.random_sinf = 0.0
        self.random_cosf = 0.0
        self.iz = 0
        self.nscat = 0
        self.chi = -1
        self.ichi = -1

    def save(self, w, nscat, iz, chi, ichi):

        self.weight = w
        self.nscat = nscat
        self.iz = iz
        self.chi = chi
        self.ichi = ichi

        return ERRCODE.SUCCESS

    def getRandomCircle(self, r2min):

        rand = Random()

        # rejection method
        for i in range(1, 201):
            w1 = 1.0 - 2.0 * float(rand.getRandom())
            w2 = 1.0 - 2.0 * float(rand.getRandom())

            self.random_r2 = w1 ** 2 + w2 ** 2

            if ((self.random_r2 > r2min) and (self.random_r2 <= 1.0)):
                self.random_r = sqrt(self.random_r2)
                self.random_cosf = w1 / self.random_r
                self.random_sinf = w2 / self.random_r
                return ERRCODE.SUCCESS

        w1 = 1.0
        w2 = 0.0
        self.random_r2 = 1.0

        # azimuth
        self.random_r = sqrt(self.random_r2)
        self.random_cosf = w1 / self.random_r
        self.random_sinf = w2 / self.random_r

        return ERRCODE.SUCCESS

    # Scattering (renewal of the direction).
    def artscat(self, vectCoord, sin_q, cos_q, sin_f, cos_f):
        MIN_VALUE = 1.0e-15
        sin_q0 = sqrt(vectCoord.x ** 2 + vectCoord.y ** 2)

        if (sin_q0 > MIN_VALUE):
            temp = 1.0 / sin_q0
            cos_f0 = vectCoord.x * temp
            sin_f0 = vectCoord.y * temp

            vectCoord.x = cos_q * vectCoord.x + sin_q * (cos_f * vectCoord.z * cos_f0 - sin_f * sin_f0)
            vectCoord.y = cos_q * vectCoord.y + sin_q * (cos_f * vectCoord.z * sin_f0 + sin_f * cos_f0)
            vectCoord.z = cos_q * vectCoord.z - sin_q * cos_f * sin_q0

        elif (vectCoord.z > 0.0):
            vectCoord.x = sin_q * cos_f
            vectCoord.y = sin_q * sin_f
            vectCoord.z = cos_q

        else:
            vectCoord.x = -sin_q * cos_f
            vectCoord.y = -sin_q * sin_f
            vectCoord.z = -cos_q

        vectCoord.nrmlVector()

        return ERRCODE.SUCCESS


    # *********************************************************************
    # Search grid number by binary search method.
    #  grd(i) <= dat < grd(i+1) or grd(i) > dat >= grd(i+1).
    #
    # Note :
    #  irvrs = 0    : grd(i) is increasing with i
    #          1    : grd(i) is decreasing with i
    # *********************************************************************
    def ifunc_vctrbinsrch(self, grd, n, irvrs, dat):

        i0 = 0
        i1 = n + 1

        if (irvrs <= 0):
            while (i1 > i0 + 1):
                i = int((i0 + i1) / 2)
                if (dat >= grd[i]):
                    i0 = i
                else:
                    i1 = i
            result = i0
        else:
            while (i1 > i0 + 1):
                i = int((i0 + i1) / 2)
                if (dat < grd[i]):
                    i0 = i
                else:
                    i1 = i
            result = i0

        return result

    # ***********************************************************************
    # Dual-end truncation approximation (DTA):
    #   Truncation of forward peak of phase function with a correction
    #   for conservativation of 1st and 2nd moments of cosQ.
    #
    # Note:
    #  1, cump() is in the reverse order (1 to 0)
    #  2, The renormalized phase function is
    #     P = phsf(ang) / (1 - ftf - ftb) if angtf < ang < angtb
    #           0                           else
    # ***********************************************************************
    def phsftrunc5(self, cosa, cump, pdf, ang, nang):

        # no truncation
        if (self.fd <= 0.0001):
            self.ftf = 0.0
            self.ftb = 0.0
            self.g1t = self.g1
            self.g2t = self.g2
            self.angtf = 0.0
            self.angtb = 180.0
            self.iangtf = 1
            self.iangtb = nang
            return ERRCODE.SUCCESS

        # initializations
        fdc = 1.0 - self.fd
        self.iangtf = self.ifunc_vctrbinsrch(cump, nang, 1, fdc)
        self.iangtb = nang
        self.g1t = (self.g1 - self.fd) / (1.0 - self.fd)
        self.g2t = (self.g2 - self.fd) / (1.0 - self.fd)
        self.g1t = max(-1.0, min(1.0, self.g1t))
        self.g2t = max(0.0, min(1.0, self.g2t))
        sum0 = sum1 = sum2 = 0.0
        g1tlo = self.g1
        g1thi = self.g1
        g2tlo = self.g2
        g2thi = self.g2
        c0 = c1 = 0.0

        # integration over the delta part
        for iAng in range(1, self.iangtf + 1):
            w = pdf[iAng]
            c0 = cosa[iAng]
            c1 = cosa[iAng + 1]
            sum0 += w
            sum1 += w * 0.5 * (c0 + c1)
            sum2 += w * 0.5 * (c0 ** 2 + c1 ** 2)
            if (iAng >= self.iangtf - 1):
                g1tlo = g1thi
                g1thi = (self.g1 - sum1) / (1.0 - sum0)

        iexit = -1

        while (iexit != 1):
            # truncation of backward region
            if (iexit == 0):
                if ((self.iangtb < nang) and (g2tlo <= self.g2t)):
                    iexit = 1
                    if (abs(g2thi - g2tlo) < 1.0e-6):
                        ratb = 0.0
                    else:
                        ratb = (self.g2t - g2tlo) / (g2thi - g2tlo)

                    self.angtb = ang[self.iangtb] * (1.0 - ratb) + ang[self.iangtb + 1] * ratb
                    self.ftb = cump[self.iangtb + 1] + pdf[self.iangtb] * (1.0 - ratb)
                    w = -pdf[self.iangtb] * ratb
                    c0 = cos(radians(self.angtb))
                    c1 = cosa[self.iangtb]
                else:
                    g2thi = g2tlo
                    self.iangtb -= 1
                    self.angtb = ang[self.iangtb]
                    self.ftb = cump[self.iangtb]
                    w = pdf[self.iangtb]
                    c0 = cosa[self.iangtb]
                    c1 = cosa[self.iangtb + 1]

                sum0 += w
                sum1 += w * 0.5 * (c0 + c1)
                sum2 += w * 0.5 * (c0 ** 2 + c1 ** 2)
                g1thi = (self.g1 - sum1) / (1.0 - sum0)

                if (sum0 > 0.99):
                    iexit = 1

                w = pdf[self.iangtf]
                c0 = cosa[self.iangtf]
                c1 = cosa[self.iangtf + 1]
                sum0a = sum0 - w
                sum1a = sum1 - w * 0.5 * (c0 + c1)
                g1tlo = (self.g1 - sum1a) / (1.0 - sum0a)
            else:
                iexit = 0

            # truncation of forward region
            while (self.g1t < g1thi):

                self.iangtf += 1
                g1tlo = g1thi
                w = pdf[self.iangtf]
                c0 = cosa[self.iangtf]
                c1 = cosa[self.iangtf + 1]
                sum0 += w
                sum1 += w * 0.5 * (c0 + c1)
                sum2 += w * 0.5 * (c0 ** 2 + c1 ** 2)
                g1thi = (self.g1 - sum1) / (1.0 - sum0)

                if (sum0 > 0.99):
                    self.angtf = ang[self.iangtf]
                    self.ftf = 1.0 - cump[self.iangtf]
                    return ERRCODE.SUCCESS

            while (self.g1t > g1tlo):
                self.iangtf -= 1
                g1thi = g1tlo
                w = -pdf[self.iangtf]
                c0 = cosa[self.iangtf]
                c1 = cosa[self.iangtf + 1]
                sum0 += w
                sum1 += w * 0.5 * (c0 + c1)
                sum2 += w * 0.5 * (c0 ** 2 + c1 ** 2)
                g1tlo = (self.g1 - sum1) / (1.0 - sum0)
                if (self.iangtf <= 1):
                    self.angtf = 0.0
                    self.ftf = 0.0
                    return ERRCODE.SUCCESS

            if (abs(g1thi - g1tlo) < 1.0e-6):
                ratf = 0.0
            else:
                ratf = (self.g1t - g1tlo) / (g1thi - g1tlo)

            self.angtf = ang[self.iangtf] * (1.0 - ratf) + ang[self.iangtf + 1] * ratf
            self.ftf = 1.0 - (cump[self.iangtf + 1] + pdf[self.iangtf] * (1.0 - ratf))
            w = pdf[self.iangtf] * (1.0 - ratf)
            c0 = cos(radians(self.angtf))
            c1 = cosa[self.iangtf + 1]

            sum0a = sum0 - w
            sum2a = sum2 - w * 0.5 * (c0 ** 2 + c1 ** 2)
            g2tlo = (self.g2 - sum2a) / (1.0 - sum0a)

        return ERRCODE.SUCCESS

    # *********************************************************************
    # Get asymmetry factors of scattering phase function for each
    #  forward and backward regions.
    #
    # Note: ang() should be an increasing function.
    # *********************************************************************
    def phsfbiasym2(self, ang, phs, nAng, angs):

        # separation point & P interpolation
        iAngs = self.ifunc_vctrbinsrch(ang, nAng, 0, angs)
        iAngs = max(1, min(nAng - 1, iAngs))
        dang = ang[iAngs + 1] - ang[iAngs]

        if (dang > 1.0e-35):
            phss = phs[iAngs]
        else:
            rat = (angs - ang[iAngs]) / dang
            phss = phs[iAngs] * (1.0 - rat) + phs[iAngs + 1] * rat

        # forward region
        sum0 = sum1 = sum2 = 0.0
        a0 = radians(ang[1])
        sin_a0 = sin(a0)
        cos_a0 = cos(a0)

        for iAng in range(1, iAngs):
            a1 = radians(ang[iAng + 1])
            sin_a1 = sin(a1)
            cos_a1 = cos(a1)
            w = (a1 - a0) * (phs[iAng + 1] * sin_a1 + phs[iAng] * sin_a0)
            sum0 += w
            sum1 += w * 0.5 * (cos_a0 + cos_a1)
            sum2 += w * 0.5 * (cos_a0 ** 2 + cos_a1 ** 2)
            a0 = a1
            sin_a0 = sin_a1
            cos_a0 = cos_a1

        a1 = radians(angs)
        sin_a1 = sin(a1)
        cos_a1 = cos(a1)
        w = (a1 - a0) * (phss * sin_a1 + phs[iAngs] * sin_a0)
        sum0 += w
        sum1 += w * 0.5 * (cos_a0 + cos_a1)
        sum2 += w * 0.5 * (cos_a0 ** 2 + cos_a1 ** 2)
        self.g0fwd = sum0
        if (sum0 > 1.0e-35):
            self.g1fwd = sum1 / sum0
            self.g2fwd = sum2 / sum0
        else:
            self.g1fwd = 0.0
            self.g2fwd = 0.0

        # backward region
        sum0 = sum1 = sum2 = 0.0
        a1 = radians(ang[nAng])
        sin_a1 = sin(a1)
        cos_a1 = cos(a1)

        for iAng in range(nAng - 1, iAngs, -1):
            a0 = radians(ang[iAng])
            sin_a0 = sin(a0)
            cos_a0 = cos(a0)
            w = (a1 - a0) * (phs[iAng + 1] * sin_a1 + phs[iAng] * sin_a0)
            sum0 += w
            sum1 += w * 0.5 * (cos_a0 + cos_a1)
            sum2 += w * 0.5 * (cos_a0 ** 2 + cos_a1 ** 2)
            a1 = a0
            sin_a1 = sin_a0
            cos_a1 = cos_a0

        a0 = radians(angs)
        sin_a0 = sin(a0)
        cos_a0 = cos(a0)
        w = (a1 - a0) * (phs[iAngs + 1] * sin_a1 + phss * sin_a0)
        sum0 += w
        sum1 += w * 0.5 * (cos_a0 + cos_a1)
        sum2 += w * 0.5 * (cos_a0 ** 2 + cos_a1 ** 2)
        self.g0bwd = sum0
        if (sum0 > 1.0e-35):
            self.g1bwd = sum1 / sum0
            self.g2bwd = sum2 / sum0
        else:
            self.g1bwd = 0.0
            self.g2bwd = 0.0

        # moments
        sumNum = self.g0fwd + self.g0bwd
        if (sumNum > 1.0e-35):
            self.g0fwd /= sumNum
            self.g0bwd /= sumNum
        else:
            self.g0fwd = 0.0
            self.g0bwd = 0.0

        self.g1 = self.g0fwd * self.g1fwd + self.g0bwd * self.g1bwd
        self.g2 = self.g0fwd * self.g2fwd + self.g0bwd * self.g2bwd

        return ERRCODE.SUCCESS

    # **********************************************************************
    # Find a grid number i, where
    #       xgrd(i) <= x  < xgrd(i+1) when xgrd is increasing vector, or
    #       xgrd(i)  > x >= xgrd(i+1) when xgrd is decreasing vector,
    #   by sequential search method with a given initial estimate, ix.
    # Output i will be in the range [ixmin, ixmax]. x should follow
    #       xgrd(ixmin) <= x  < xgrd(ixmax+1), or
    #       xgrd(ixmin)  > x >= xgrd(ixmax+1).
    # **********************************************************************
    def i_rvctrssrch(self, xGrd, x, ix, ixMin, ixMax):

        result = -1

        # Increasing vector
        if ((xGrd[ixMin]) < (xGrd[ixMax + 1])):
            if (x < xGrd[ix]):
                for i in range(ix - 1, ixMin - 1, -1):
                    if (x >= xGrd[i]):
                        result = i
                        break
                result = max(ixMin, result)
            else:
                for i in range(ix + 1, ixMax + 2):
                    if (x < xGrd[i]):
                        result = i
                        break
                result = min(ixMax, result - 1)
        # Decreasing vector
        else:
            if (x >= xGrd[ix]):
                for i in range(ix - 1, ixMin - 1, -1):
                    if (x < xGrd[i]):
                        result = i
                        break
                result = max(ixMin, result)
            else:
                for i in range(ix + 1, ixMax + 2):
                    if (x >= xGrd[i]):
                        result = i
                        break
                result = min(ixMax, result - 1)

        return result

    # ************************************************************************
    # Make LUTs of scattering angle for uniform distributions and others.
    #
    # Notes
    # wrkang: interpolated angles (degree)
    # wrkphs: interpolated phase functions
    # wrkpdf: normalized PDFs for spherical integration
    # wrkcum: normalized cumulative PDFs by backward integration (pi to 0)
    # wrksca: scattering angle LUT for uniform distribution
    # ************************************************************************
    def artsanglut(self, wrkang, wrkphs, wrkcum, wrkpdf, wrksca, nwrk, nlut):

        # Integrate P (backward)
        sumNum = 0.0
        wrkpdf[nwrk] = 0.0
        wrkcum[nwrk] = 0.0
        a1 = radians(wrkang[nwrk])
        sin_a1 = sin(a1)

        for iwrk in range(nwrk - 1, 0, -1):
            a0 = radians(wrkang[iwrk])
            sin_a0 = sin(a0)
            p = (a1 - a0) * (wrkphs[iwrk + 1] * sin_a1 + wrkphs[iwrk] * sin_a0)
            sumNum += p

            wrkpdf[iwrk] = p
            wrkcum[iwrk] = sumNum
            a1 = a0
            sin_a1 = sin_a0

        # Normalizations
        aSum = 1.0 / sumNum
        f = 4.0 / sumNum
        for iwrk in range(1, nwrk + 1):
            wrkphs[iwrk] *= f
            wrkcum[iwrk] *= aSum
            wrkpdf[iwrk] *= aSum
        wrkcum[1] = 1.0
        wrkcum[nwrk] = 0.0

        # scattering angle
        pdlt = 1.0 / float(nlut)
        iwrk = 1
        for iLut in range(1, nlut):
            p = pdlt * float(nlut - iLut)
            iwrk = self.i_rvctrssrch(wrkcum, p, iwrk, iwrk, nwrk - 1)
            if (wrkpdf[iwrk] > 1.0e-6):
                rat = (p - wrkcum[iwrk + 1]) / wrkpdf[iwrk]
                wrksca[iLut] = rat * wrkang[iwrk] + (1.0 - rat) * wrkang[iwrk + 1]
            else:
                wrksca[iLut] = wrkang[iwrk]

        wrksca[0] = 0.0
        wrksca[nlut] = 180.0

        return ERRCODE.SUCCESS

    # ########################################################
    # Interpolation of phase function by natural cubic spline
    # rawang is angle, rawphs is phs function
    # wrkc1, wrkc2, wrkc3 are empty
    # wrkphs is empty
    # nraw is the number of angle. user input
    # nwrk is the number of LUT
    # ########################################################
    def artphsfintp(self, rawang, rawphs, wrkc1, wrkc2, wrkc3, wrkang, wrkphs, nraw, nwrk):

        wrkc1[1] = (rawphs[2] - rawphs[1]) / (rawang[2] - rawang[1])
        wrkc1[2] = (rawphs[3] - rawphs[2]) / (rawang[3] - rawang[2])
        wrkc2[2] = 2.0 * (rawang[3] - rawang[1])
        wrkc3[2] = wrkc1[2] - wrkc1[1]

        for i in range(3, nraw):
            h1 = rawang[i + 1] - rawang[i]
            h0 = rawang[i] - rawang[i - 1]

            temp = h0 / wrkc2[i - 1]
            wrkc1[i] = (rawphs[i + 1] - rawphs[i]) / h1
            wrkc2[i] = 2.0 * (h0 + h1) - h0 * temp
            wrkc3[i] = wrkc1[i] - wrkc1[i - 1] - wrkc3[i - 1] * temp

        sig0 = wrkc3[nraw - 1] / wrkc2[nraw - 1]
        h = rawang[nraw] - rawang[nraw - 1]
        wrkc1[nraw - 1] -= h * 2.0 * sig0
        wrkc2[nraw - 1] = 3.0 * sig0
        wrkc3[nraw - 1] = -sig0 / h
        sig1 = sig0

        for i in range(nraw - 2, 1, -1):
            h = rawang[i + 1] - rawang[i]
            sig0 = (wrkc3[i] - h * sig1) / wrkc2[i]
            wrkc1[i] -= h * (sig1 + 2.0 * sig0)
            wrkc2[i] = 3.0 * sig0
            wrkc3[i] = (sig1 - sig0) / h
            sig1 = sig0

        h1 = rawang[2] - rawang[1]
        wrkc1[1] -= h1 * sig1
        wrkc2[1] = 0.0
        wrkc3[1] = sig1 / h1

        # Interpolation
        iraw = 1
        for iwrk in range(1, nwrk + 1):
            x = wrkang[iwrk]
            iraw = self.i_rvctrssrch(rawang, x, iraw, iraw, nraw - 1)
            dx = x - rawang[iraw]
            y = rawphs[iraw] + (wrkc1[iraw] + (wrkc2[iraw] + wrkc3[iraw] * dx) * dx) * dx
            wrkphs[iwrk] = max(0.0, y)

        return ERRCODE.SUCCESS

    # Normalize phase function to -sumnorm-.
    def nrmlzphsf(self, nAng, ang, pf, sumNorm):

        sumNum = 0.0

        for i in range(nAng - 1, 0, -1):
            a0 = radians(ang[i])
            a1 = radians(ang[i + 1])
            sumNum += (a0 - a1) * 0.5 * (pf[i + 1] * sin(a1) + pf[i] * sin(a0))

        sumNum = abs(sumNum)

        f = 2.0 * sumNorm / sumNum

        for i in range(1, nAng + 1):
            pf[i] *= f

        return ERRCODE.SUCCESS

    # Horizontal shift
    def arthshift(self, x, ux, path, xmax):

        x += path * ux

        if ((x < 0.0) or (x >= xmax)):
            x = fmod(x, xmax)
            # logging.debug("In arthshift, x <0, x = " + str(x) + ", before x = " + str(t))
            if (x < 0.0):
                x += xmax

        return x


    # Make LUTs for phase function & scattering angle.
    def optics(self, ext, omg, phs, angle, nmix):
        kNRaw = 5000
        nLut1 = comm.N_LUT1 + 1
        nLut = comm.N_LUT

        rawAngle = [0.0] * kNRaw
        rawPhs = [0.0] * kNRaw
        wrkC1 = [0.0] * kNRaw
        wrkC2 = [0.0] * kNRaw
        wrkC3 = [0.0] * kNRaw

        wrkAngle = [0.0] * nLut1
        wrkPhs = [0.0] * nLut1
        wrkCum = [0.0] * nLut1
        wrkPdf = [0.0] * nLut1
        wrkCos = [0.0] * nLut1
        wrkSca = [0.0] * nLut1

        # Initializations
        comm.FS_ANG = float(nLut - 1) / pi
        nChi = 6
        chiHi = 0.9
        chiLo = 0.4
        fMax = 0.8
        chiBin = (chiHi - chiLo) / float(nChi - 2)
        for iChi in range(1, nChi):
            comm.CHI_GRD[iChi] = chiHi - chiBin * float(iChi - 1)
        comm.CHI_GRD[nChi] = -1.0

        # angle & cosine
        delt = 180.0 / float(nLut)
        nwrk = nLut
        # print()
        for iWrk in range(1, nwrk + 1):
            temp = delt * (iWrk - 1)
            wrkAngle[iWrk] = temp
            wrkCos[iWrk] = cos(radians(temp))

        # Normalization of the phase functions
        for iMix in range(1, nmix + 1):
            for iAng in range(1, comm.N_ANG + 1):
                rawPhs[iAng] = phs[iMix, iAng]

            self.nrmlzphsf(comm.N_ANG, angle, rawPhs, 1.0)

            for iAng in range(1, comm.N_ANG + 1):
                phs[iMix, iAng] = rawPhs[iAng]

        # loop for layers
        for iz in range(1, comm.N_Z + 1):

            # mix optical properties
            sumA = 0.0
            sumE = 0.0
            sumS = 0.0

            rawPhs = [0.0] * (comm.N_ANG + 1)

            for iMix in range(1, nmix + 1):
                ee = ext[iMix, iz]
                o = omg[iMix]
                s = ee * o
                sumA += ee * (1.0 - o)
                sumE += ee
                sumS += s

                for iAng in range(1, comm.N_ANG + 1):
                    rawPhs[iAng] += s * phs[iMix, iAng]

            comm.EXT_T1D[iz, 1] = sumE
            comm.ABS_T1D[iz] = sumA

            f = 1.0 / max(1.0e-35, sumS)

            for iAng in range(1, comm.N_ANG + 1):
                rawAngle[iAng] = angle[iAng]
                rawPhs[iAng] *= f

            nRaw = comm.N_ANG

            # make LUT
            self.artphsfintp(rawAngle, rawPhs, wrkC1, wrkC2, wrkC3, wrkAngle, wrkPhs, nRaw, nwrk)
            self.artsanglut(wrkAngle, wrkPhs, wrkCum, wrkPdf, wrkSca, nwrk, nLut - 1)

            for iLut in range(nLut):
                #print(iLut)
                comm.PF_LUT[iLut, iz] = wrkPhs[iLut + 1]
                comm.SA_LUT[iLut, iz] = radians(wrkSca[iLut])

            comm.PF_LUT[comm.N_LUT1 - 1, iz] = comm.PF_LUT[nLut - 1, iz]
            comm.SA_LUT[comm.N_LUT1 - 1, iz] = comm.PF_LUT[nLut - 1, iz]

            # truncation
            self.phsfbiasym2(wrkAngle, wrkPhs, nwrk, 90.0)

            cftMAX = fMax * self.g0fwd * self.g1fwd ** 4 / float(nChi - 1)

            ee = comm.EXT_T1D[iz, 1]
            s = ee - comm.ABS_T1D[iz]

            for iChi in range(1, nChi + 1):
                self.fd = cftMAX * float(iChi - 1)

                if (self.fd <= 0.0):
                    comm.TRU_LUT[1, iChi, iz] = 1.0
                    comm.TRU_LUT[2, iChi, iz] = 0.0
                    comm.TRU_LUT[3, iChi, iz] = 1.0
                    comm.TRU_LUT[4, iChi, iz] = -1.0e35
                    comm.TRU_LUT[5, iChi, iz] = 1.0e35
                    comm.TRU_LUT[6, iChi, iz] = abs(self.g1)

                else:
                    self.phsftrunc5(wrkCos, wrkCum, wrkPdf, wrkAngle, nwrk)

                    ftt = 1.0 - (self.ftf + self.ftb)
                    comm.TRU_LUT[1, iChi, iz] = 1.0 - self.fd
                    comm.TRU_LUT[2, iChi, iz] = self.ftf
                    comm.TRU_LUT[3, iChi, iz] = ftt
                    comm.TRU_LUT[4, iChi, iz] = radians(self.angtf)
                    comm.TRU_LUT[5, iChi, iz] = radians(self.angtb)
                    comm.TRU_LUT[6, iChi, iz] = abs(self.g1t)

                comm.EXT_T1D[iz, iChi] = ee - s * self.fd

        return ERRCODE.SUCCESS

    # Traces a trajectory in plane-parallel vertically-inhomogeneous atmosphere
    def escape(self, phoCoord, vectCoord, iz, ichi, ikd, resultCoord):

        resultCoord.setPosition(phoCoord.x, phoCoord.y, phoCoord.z)
        # logging.debug("MC1D escape start...")
        # TAU integration
        tau = 0.0
        if (ikd == 0):
            logging.debug("MC1D escape finish (ikd = 0).")
            return ERRCODE.SUCCESS

        if (vectCoord.z >= 0.0):
            for izr in range(iz, comm.N_Z + 1):
                tau += (comm.Z_GRD[izr] - resultCoord.z) * (comm.ABS_G1D[izr, ikd] + comm.EXT_T1D[izr, ichi])
                resultCoord.z = comm.Z_GRD[izr]

            tau /= vectCoord.z
        else:
            for izr in range(iz, 0, -1):
                tau += (resultCoord.z - comm.Z_GRD[izr - 1]) * (comm.ABS_G1D[izr, ikd] + comm.EXT_T1D[izr, ichi])
                resultCoord.z = comm.Z_GRD[izr - 1]

            tau /= (-1.0 * vectCoord.z)

        # escape location
        fpath = (resultCoord.z - phoCoord.z) / abs(vectCoord.z)
        resultCoord.x = self.arthshift(resultCoord.x, vectCoord.x, fpath, comm.X_MAX)
        resultCoord.y = self.arthshift(resultCoord.y, vectCoord.y, fpath, comm.Y_MAX)
        # logging.debug("MC1D escape finish.")
        return tau

    # Traces a trajectory in plane-parallel vertically-inhomogeneous atmosphere
    def trace(self, phoCoord, vectCoord, w, wq, ftau, chi, ikd, iz, nscat, ichi):

        rand = Random()

        logging.debug("Monte Carlo 1-D simulation start...")
        string = "x = " + str(phoCoord.x) + ", y =" + str(phoCoord.y) + ", z =" + str(phoCoord.z)
        logging.debug(string)
        # loop for layers
        while ((iz >= 1) and (iz <= comm.N_Z)):
            string = "x = " + str(phoCoord.x) + ", y =" + str(phoCoord.y) + ", z =" + str(phoCoord.z)
            logging.debug("mc1dtrace: start while")
            logging.debug(string)
            logging.debug("uz = " + str(vectCoord.z))
            absg = comm.ABS_G1D[iz, ikd]
            extm = absg + comm.EXT_T1D[iz, ichi]
            absm = absg + comm.ABS_T1D[iz]

            if (vectCoord.z < 0.0):
                path = (comm.Z_GRD[iz - 1] - phoCoord.z) / vectCoord.z
            else:
                path = (comm.Z_GRD[iz] - phoCoord.z) / vectCoord.z
            logging.debug("vectCoord.z = " + str(vectCoord.z))
            logging.debug("path = " + str(path))
            # scattering events
            if (extm > 1.0e-20):
                fpath = ftau / extm

                while(fpath <= path):
                    # motion to the collision point
                    phoCoord.z += fpath * vectCoord.z

                    phoCoord.x = self.arthshift(phoCoord.x, vectCoord.x, fpath, comm.X_MAX)
                    phoCoord.y = self.arthshift(phoCoord.y, vectCoord.y, fpath, comm.Y_MAX)
                    string = "x = " + str(phoCoord.x) + ", y =" + str(phoCoord.y) + ", z =" + str(phoCoord.z)
                    logging.debug("mc1dtrace: after inner while")
                    logging.debug(string)
                    logging.debug("uz = " + str(vectCoord.z))
                    # weight scaling
                    w *= (1.0 - absm / extm)

                    # Russian roulette: survive or killed?
                    if (w < comm.WRR * 0.5):
                        if (w > comm.WRR * rand.getRandom()):
                            w = comm.WRR    # survive
                        else:
                            w = 0.0         # killed

                    if (w <= 0.0):
                        self.save(w, nscat, iz, chi, ichi)
                        logging.debug("Monte Carlo 1-D simulation finish.(w <= 0)")
                        return ERRCODE.LOW_WEIGHT

                    # properties
                    ftf = comm.TRU_LUT[2, ichi, iz]
                    ftt = comm.TRU_LUT[3, ichi, iz]
                    utf = comm.TRU_LUT[4, ichi, iz]
                    utb = comm.TRU_LUT[5, ichi, iz]
                    chi *= comm.TRU_LUT[6, ichi, iz]

                    nscat += 1
                    # print("ichi:", ichi)
                    while (chi < comm.CHI_GRD[ichi]):
                        ichi += 1

                    # local estimates
                    for irdc in range(1, comm.N_RDC + 1):
                        vectR = Position()
                        vectR.setPosition(comm.UX_RTAB[irdc],
                                          comm.UY_RTAB[irdc],
                                          comm.UZ_RTAB[irdc])

                        vectR.angle_twoVectors(vectCoord)
                        q = vectR.angle
                        cos_q = vectR.cosA

                        if ((q < utf) or (q > utb)):
                            continue

                        rilut = comm.FS_ANG * q
                        ilut = int(rilut)
                        rat = rilut - float(ilut)
                        pf = (1.0 - rat) * comm.PF_LUT[ilut, iz] + rat * comm.PF_LUT[ilut + 1, iz]
                        adf = pf / (4.0 * pi * abs(vectR.z) * ftt)

                        resultCoord = Position()
                        tau = self.escape(phoCoord, vectR, iz, ichi, ikd, resultCoord)

                        # give here the # of pixels
                        ixr = int(resultCoord.x / comm.X_MAX * comm.K_NXR + 1)
                        iyr = int(resultCoord.y / comm.Y_MAX * comm.K_NYR + 1)

                        p = w * adf * exp(-tau)

                        ixr = min(ixr, comm.K_NXR)
                        iyr = min(iyr, comm.K_NYR)

                        comm.PROC_F[ixr, iyr, irdc] += p
                        comm.PROC_Q[ixr, iyr, irdc] += p * wq

                    for i in range(1, comm.N_ANG_C + 1):
                        vectR = Position()
                        vectR.setPosition(comm.URC_coord[i].x,
                                          comm.URC_coord[i].y,
                                          comm.URC_coord[i].z)

                        vectR.angle_twoVectors(vectCoord)
                        q = vectR.angle
                        cos_q = vectR.cosA

                        if ((q < utf) or (q > utb)):
                            continue

                        rilut = comm.FS_ANG * q
                        ilut = int(rilut)
                        rat = rilut - float(ilut)
                        pf = (1.0 - rat) * comm.PF_LUT[ilut, iz] + rat * comm.PF_LUT[ilut + 1, iz]
                        adf = pf / (4.0 * pi * abs(vectR.z) * ftt)

                        resultCoord = Position()
                        tau = self.escape(phoCoord, vectR, iz, ichi, ikd, resultCoord)

                        ixr = resultCoord.x / comm.X_MAX * comm.K_NXR + 1
                        iyr = resultCoord.y / comm.Y_MAX * comm.K_NYR + 1

                        p = w * adf * exp(-tau)

                        # brf(2, i) = 0 in the previous stage
                        comm.BRF[2, i] += p

                    # scattering
                    self.getRandomCircle(1.0e-12)

                    rlut = (comm.N_LUT - 1) * (ftf + ftt * self.random_r2)
                    ilut = int(rlut)
                    rat = rlut - float(ilut)

                    s = (1.0 - rat) * comm.SA_LUT[ilut, iz] + rat * comm.SA_LUT[ilut + 1, iz]
                    cos_q = cos(s)
                    sin_q = sin(s)

                    self.artscat(vectCoord, sin_q, cos_q, self.random_sinf, self.random_cosf)

                    if (abs(vectCoord.z) < 1.0e-17):
                        vectCoord.z = 1.0e-17

                    # new path
                    ftau = - log(max(1.0e-35, float(rand.getRandom())))

                    if (vectCoord.z < 0.0):
                        path = (comm.Z_GRD[iz - 1] - phoCoord.z) / vectCoord.z
                    else:
                        path = (comm.Z_GRD[iz] - phoCoord.z) / vectCoord.z

                    extm = absg + comm.EXT_T1D[iz, ichi]
                    fpath = ftau / extm

                ftau = max(0.0, ftau - extm * path)

            # bounds
            if (vectCoord.z < 0.0):
                # downward
                iz -= 1
                phoCoord.z = comm.Z_GRD[iz]
            else:
                # upward
                phoCoord.z = comm.Z_GRD[iz]
                iz += 1

            phoCoord.x = self.arthshift(phoCoord.x, vectCoord.x, path, comm.X_MAX)
            phoCoord.y = self.arthshift(phoCoord.y, vectCoord.y, path, comm.Y_MAX)
            string = "x = " + str(phoCoord.x) + ", y =" + str(phoCoord.y) + ", z =" + str(phoCoord.z)
            logging.debug("mc1dtrace: after while")
            logging.debug(string)
            logging.debug("uz = " + str(vectCoord.z))

        self.save(w, nscat, iz, chi, ichi)
        logging.debug("Monte Carlo 1-D simulation finish.")
        return ERRCODE.SUCCESS


