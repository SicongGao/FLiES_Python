import random
from math import *
from iparam import Parameters


class MonteCarlo:

    nscat = 0
    cNscat = 0
    cIchi = 0
    cIkd = 0

    def __init__(self):
        return

    def stem(self, w, wq, phoCoord, vectCoord, nscat, tObj, face, ichi, ikd, iParameter):
        return

    def canopy(self):
        return

    def floor(self):
        return







# ###################
# FOR TEST
# ###################


def test():
    R = 0.5
    inside = 0
    allNum = 10000
    for i in range(allNum):
        rX = random.uniform(0, 1)
        rY = random.uniform(0, 1)
        rX = (rX - R) * (rX - R)
        rY = (rY - R) * (rY - R)
        rNEW = math.sqrt(rX + rY)

        if (rNEW <= R):
            inside += 1

    PI = 4 * inside / allNum
    print("PI = ", PI)
    print("ERR = ", abs((math.pi - PI) / math.pi))

