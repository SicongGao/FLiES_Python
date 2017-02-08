import random
import math

class MonteCarlo:

    def __init__(self):
        return

    def stem(self):
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

