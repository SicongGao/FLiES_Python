import ERRCODE
import math
import common


class G_Function:

    GT_BLC = [0.0] * 180
    GT_BLB = [0.0] * 180
    GT_BLF = [0.0] * 180
    __x = []
    __w = []

    def __init__(self):
        self.__w = [0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443]
        self.__x = [0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285]

    def igtbl(self, para):
        max = math.pi
        min = 0.0

        xm = 0.5 * (max + min)
        xr = 0.5 * (max - min)

        for i in range(180):
            th = math.radians(i)
            sc = 0.0
            sb = 0.0
            sf = 0.0

            for j in range(len(self.__w)):
                dx = xr * self.__x[j]
                # for canopy
                sc += self.__w[j] * (self.fgl(xm + dx, common.M_C) * self.fpsi(th, xm + dx) +
                                     self.fgl(xm - dx, common.M_C) * self.fpsi(th, xm - dx))
                # for branch area
                sb += self.__w[j] * (self.fgl(xm + dx, common.M_B) * self.fpsi(th, xm + dx) +
                                     self.fgl(xm - dx, common.M_B) * self.fpsi(th, xm - dx))
                # for forest floor
                sf += self.__w[j] * (self.fgl(xm + dx, common.M_F) * self.fpsi(th, xm + dx) +
                                     self.fgl(xm - dx, common.M_F) * self.fpsi(th, xm - dx))

            self.GT_BLC[i] = sc * xr
            self.GT_BLB[i] = sb * xr
            self.GT_BLF[i] = sf * xr

        return ERRCODE.SUCCESS

    #     ********************************
    #     psi function defined in Shultis and Myneni (1989)
    #     this function has been validated
    #     2008/01/16
    #     ********************************
    def fpsi(self, th, thl):
        result = 0.0
        pht = -common.getCos(th) * common.getCos(thl)

        if (() <= 1E-5):
            pht /= 1E-5
        else:
            pht /= common.getSin(th) * common.getSin(thl)

        absa = abs(pht)

        if(absa > 1.0):
            result = abs(common.getCos(th) * common.getCos(thl))
        else:
            pht = common.getACos(pht)
            result = common.getCos(th) * common.getCos(thl) * (2.0 * pht / math.pi - 1.0) + \
                     (2.0 / math.pi) * common.getSin(th) * common.getSin(thl) * common.getSin(pht)
        return ERRCODE.SUCCESS

    #     ********************************
    #     leaf angle distribution function
    #     (gl x sin(th))
    #     using Bunnik (1978) definition
    #     this function has been validated
    #     2008/01/16
    #     ********************************
    def fgl(self, thl, i):
        ag = [0.0, 1.0, 1.0]
        bg = [0.0, 1.0, -1.0]
        cg = [0.0, 1.0, 1.0]
        dg = [1.0, 0.0, 0.0]

        # Preparation of the coeffifient
        # i = 1: uniform, i = 2 planophile i = 3 erectrophile
        # I think when i = 1, it doesn't mean uniform, based on the paper,
        # it should be Spherical. Check Ross paper(Page 252).
        return (2.0 / math.pi) * (ag[i] + bg[i] * common.getCos(2.0 * cg[i] * thl)) + \
               dg[i] * common.getSin(thl)

