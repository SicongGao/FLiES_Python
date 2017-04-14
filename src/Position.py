from math import *

class Position:
    x = 0.0
    y = 0.0
    z = 0.0

    def __init__(self):
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.cosA = 0.0
        self.angle = 0.0

    def setPosition(self, x1, y1, z1):
        self.x = x1
        self.y = y1
        self.z = z1

    def clearPotion(self):
        self.__init__()

    def movePositionDistance(self, distance, vectCoord, x_max, y_max):
        mgn = 1.0e-2
        self.x += (distance + mgn) * vectCoord.x
        self.y += (distance + mgn) * vectCoord.y
        self.z += (distance + mgn) * vectCoord.z

        self.x -= (trunc(self.x / x_max) - 0.5 + copysign(0.5, self.x)) * x_max
        self.y -= (trunc(self.y / y_max) - 0.5 + copysign(0.5, self.y)) * y_max

    def movePosition(self, vectCoord, x_max, y_max):
        mgn = 1.0e-2
        self.x += mgn * vectCoord.x
        self.y += mgn * vectCoord.y
        self.z += mgn * vectCoord.z

        self.x -= (trunc(self.x / x_max) - 0.5 + copysign(0.5, self.x)) * x_max
        self.y -= (trunc(self.y / y_max) - 0.5 + copysign(0.5, self.y)) * y_max

    def moveVector(self, distance, vectCoord):
        self.x += distance * vectCoord.x
        self.y += distance * vectCoord.y
        self.z += distance * vectCoord.z

    def angle_twoVectors(self, vectCoord):
        self.cosA = self.x * vectCoord.x + self.y * vectCoord.y + self.z * vectCoord.z

        if ((self.cosA > -0.99) and (self.cosA < 0.99)):
            self.angle = acos(self.cosA)
        else:
            if (self.cosA < 0.0):
                dux = self.x + vectCoord.x
                duy = self.y + vectCoord.y
                duz = self.z + vectCoord.z
            else:
                dux = self.x - vectCoord.x
                duy = self.y - vectCoord.y
                duz = self.z - vectCoord.z

            aa = dux ** 2 + duy ** 2 + duz ** 2
            sinAA = aa * (1.0 - 0.25 * aa)
            self.angle = sqrt(sinAA) * ((6.0 - 2.0 * sinAA) / (6.0 - 3.0 - sinAA))

            if (self.cosA < 0.0):
                self.angle = pi - self.angle

    # ******************************************************************
    # Renormalize a unit vector by using Newton method.
    # A given unit vector must satisfy that
    # ux*ux + uy*uy + uz*uz is nearly equal to 1.
    # ******************************************************************
    def nrmlVector(self):
        # This is approximation by 1st-order Newton method
        fone = 0.5 * (3.0 - (self.x ** 2 + self.y ** 2 + self.z ** 2))

        self.x *= fone
        self.y *= fone
        self.z *= fone





