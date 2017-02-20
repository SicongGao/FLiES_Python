import ERRCODE
from Random import Random
from math import *
import common as comm
from iparam import Parameters
from Position import Position
import RussianRoulette
import common as comm

class MonteCarlo:

    nscat = 0
    cNscat = 0
    cIchi = 0
    cIkd = 0
    weight = -1

    def __init__(self):

        self.nscat = 0
        self.cNscat = 0
        self.cIchi = 0
        self.cIkd = 0
        self.weight = -1

        return ERRCODE.SUCCESS

    def save(self, nscat, w):
        self.cNscat = nscat
        self.weight = w

    def transformCoordinate(self, oriVecCoord, angleA, angleB, destVecCoord):

        MIN_VALUE = 1.0e-15

        sinA = sin(angleA)
        sinB = sin(angleB)
        cosA = cos(angleA)
        cosB = cos(angleB)
        sinT = sqrt(oriVecCoord.x ** 2 + oriVecCoord.y ** 2)
        sinP = oriVecCoord.y / sinT
        cosP = oriVecCoord.x / sinT

        if (oriVecCoord.x ** 2 + oriVecCoord.y ** 2 > MIN_VALUE):
            destVecCoord.x = cosA * oriVecCoord.x + sinA * (cosB * oriVecCoord.z * cosP - sinB * sinP)
            destVecCoord.y = cosA * oriVecCoord.y + sinA * (cosB * oriVecCoord.z * sinP + sinB * cosP)
            destVecCoord.z = cosA * oriVecCoord.z - sinA * cosB * sinT

        else:
            destVecCoord.x = sinA * cosB * copysign(1.0, oriVecCoord.z)
            destVecCoord.x = sinA * sinB * copysign(1.0, oriVecCoord.z)
            destVecCoord.z = cosA * copysign(1.0, oriVecCoord.z)

        # convert to unit vector
        c = destVecCoord.x ** 2 + destVecCoord.y ** 2 + destVecCoord.z ** 2
        c = sqrt(c)
        destVecCoord.x /= c
        destVecCoord.y /= c
        destVecCoord.z /= c

        return ERRCODE.SUCCESS

    def stem(self, w, wq, phoCoord, vectCoord, nscat, tObj, face, str, ichi, ikd, iParameter):
        MGN = 1.0e-2
        MIN_VALUE = 1.0e-8
        FD = 0.0
        CB = 3
        UZM = 0.0174524
        a = 0

        randomMethod = Random()

        # reflectance at the side of stem
        if (face == 1):
            # stem normal vector
            rx = phoCoord.x - tObj[1]
            ry = phoCoord.y - tObj[2]
            rr = sqrt(rx ** 2 + ry ** 2)

            tempCoord = Position()
            tempCoord.setPosition(rx / rr, ry / rr, 0)

            th = 0.5 * acos(1.0 - 2.0 * randomMethod.getRandom())
            ph = 2.0 * pi * randomMethod.getRandom()

            a = acos(tempCoord.x)
            cb = 6

            self.transformCoordinate(tempCoord, th, ph, vectCoord)

            if (abs(vectCoord.z) < UZM):
                vectCoord.z = copysign(UZM, vectCoord.z)

        # reflectance at the bottom of stem
        elif (face == 2):
            th = 0.5 * pi + acos(1.0 - 2.0 * randomMethod.getRandom())
            ph = 2.0 * pi * randomMethod.getRandom()
            vectCoord.setPosition(sin(th) * cos(ph),
                                  sin(th) * sin(ph),
                                  cos(th))
            phoCoord.z -= MGN

            if (abs(vectCoord.z) < UZM):
                vectCoord.z = copysign(UZM, vectCoord.z)

            a = 1.0

        # reflectance at the top of stem
        else:
            th = 0.5 * acos(1.0 - 2.0 * randomMethod.getRandom())
            ph = 2.0 * pi * randomMethod.getRandom()
            vectCoord.setPosition(sin(th) * cos(ph),
                                  sin(th) * sin(ph),
                                  cos(th))
            phoCoord.z += MGN

            if (abs(vectCoord.z) < UZM):
                vectCoord.z = copysign(UZM, vectCoord.z)

            a = 1.0

        # fpar samping (leave or branch)
        # every tree specie may have the same reflectance (str)
        comm.B_FPR += w * wq * (1.0 - str)
        # ipara.AP_NP range (0,100)
        iz = int(phoCoord.z) + 1
        if (iz > 99):
            print("Z overhight: MonteCarlo.py")
            return ERRCODE.OUT_OF_RANGE
        else:
            iParameter.AP_NP[iz] += w * wq * (1.0 - str)

        w *= str
        nscat += 1

        # Russian Roulette
        w = RussianRoulette.roulette(w)

        # call vegrad()



        self.save(nscat, w)
        return ERRCODE.SUCCESS

    def canopy(self):
        return

    def floor(self):
        return







# ##############################
# FOR TEST
# Monte Carlo in Calculating PI
# ##############################
#
# def test():
#     R = 0.5
#     inside = 0
#     allNum = 10000
#     for i in range(allNum):
#         rX = random.uniform(0, 1)
#         rY = random.uniform(0, 1)
#         rX = (rX - R) * (rX - R)
#         rY = (rY - R) * (rY - R)
#         rNEW = math.sqrt(rX + rY)
#
#         if (rNEW <= R):
#             inside += 1
#
#     PI = 4 * inside / allNum
#     print("PI = ", PI)
#     print("ERR = ", abs((math.pi - PI) / math.pi))
#
# print(3 ** 2)